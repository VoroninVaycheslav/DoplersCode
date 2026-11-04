"""Independent ENDF parsing, exact polynomial positivity, and Monte Carlo checks.

Uses only the Python standard library; the collision implementation is Fortran.
"""
import bisect
import csv
import hashlib
import io
import json
import math
import pathlib
import subprocess
from fractions import Fraction

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / 'data/u238'


def number(field):
    text = field.strip()
    for i in range(len(text)-1, 0, -1):
        if text[i] in '+-':
            return float(text[:i] + 'e' + text[i:])
    return float(text or 0)


def records(name):
    return [[number(line[i:i+11]) for i in range(0,66,11)]
            for line in (DATA / name).read_text().splitlines()]


def xs(name):
    rows = records(name)
    assert int(rows[1][4]) == 1 and rows[2][:2] == [rows[1][5],2]
    flat = [v for row in rows[3:] for v in row][:2*int(rows[1][5])]
    return flat[::2], flat[1::2]


def interpolate(grid, values, energy):
    i = min(len(grid)-2, bisect.bisect_right(grid,energy)-1)
    t = (energy-grid[i])/(grid[i+1]-grid[i])
    return (1-t)*values[i] + t*values[i+1]


def angular_nodes():
    rows = records('angular.endf')
    assert rows[0][2:4] == [0,1] and rows[1][2:4] == [0,2]
    assert rows[3][:2] == [rows[2][5],2]
    nodes = []
    i = 4
    while i < len(rows):
        energy, order = rows[i][1], int(rows[i][4])
        i += 1
        count = (order+5)//6
        coeff = [v for row in rows[i:i+count] for v in row][:order]
        nodes.append((energy,coeff))
        i += count
    assert len(nodes) == int(rows[2][5])
    return nodes


def coefficients(nodes, energy):
    grid = [row[0] for row in nodes]
    order = max(len(row[1]) for row in nodes)
    padded = [c+[0.0]*(order-len(c)) for _,c in nodes]
    return [interpolate(grid,[row[l] for row in padded],energy) for l in range(order)]


def primitive(coeff, x):
    """Independent polynomial integration in ordinary powers of x."""
    power = legendre_power(coeff)
    return sum(float(c)*(x**(k+1)-(-1)**(k+1))/(k+1) for k,c in enumerate(power))


def legendre_power(coeff):
    n = len(coeff)
    power = [Fraction(1,2)] + [Fraction(0)]*n
    previous = [Fraction(1)]
    current = [Fraction(0),Fraction(1)]
    for l,a in enumerate(coeff,1):
        for k,c in enumerate(current):
            power[k] += Fraction(a)*Fraction(2*l+1,2)*c
        following = [Fraction(0)]*(l+2)
        for k,c in enumerate(current):
            following[k+1] += Fraction(2*l+1,l+1)*c
        for k,c in enumerate(previous):
            following[k] -= Fraction(l,l+1)*c
        previous,current = current,following
    return power


def certify_positive(coeff):
    """Exact rational Bernstein subdivision: certify p(mu)>=0 everywhere."""
    power = legendre_power(coeff)
    n = len(power)-1
    # Substitute mu=2*t-1, then change power basis to Bernstein basis.
    transformed = [sum(power[k]*math.comb(k,j)*2**j*(-1)**(k-j)
                       for k in range(j,n+1)) for j in range(n+1)]
    bernstein = [sum(transformed[j]*Fraction(math.comb(k,j),math.comb(n,j))
                     for j in range(k+1)) for k in range(n+1)]
    stack = [(bernstein,0)]
    while stack:
        b,depth = stack.pop()
        if min(b) >= 0:
            continue
        assert max(b) >= 0 and depth < 30, 'Cannot certify nonnegative angular PDF'
        left,right = [b[0]],[b[-1]]
        while len(b) > 1:
            b = [(x+y)/2 for x,y in zip(b,b[1:])]
            left.append(b[0])
            right.append(b[-1])
        stack.extend([(left,depth+1),(right[::-1],depth+1)])


def main():
    manifest = json.loads((DATA/'provenance.json').read_text())
    for name,digest in manifest['section_sha256'].items():
        assert hashlib.sha256((DATA/name).read_bytes()).hexdigest() == digest, name
    nodes = angular_nodes()
    # Convex energy interpolation preserves positivity between certified nodes.
    for energy,coeff in nodes:
        if 150000 <= energy <= 19000000:
            certify_positive(coeff)
    print('PASS: source-section hashes and exact-rational PDF positivity on 150 keV–19 MeV nodes')
    total_grid,total = xs('total.endf')
    elastic_grid,elastic = xs('elastic.endf')
    for energy in (1100000,18000000):
        n = 16000
        output = subprocess.check_output([str(ROOT/'build/scatter_demo'),str(energy),str(n),'76543'],
                                         cwd=ROOT,text=True)
        metadata = dict(line[2:].split(' = ') for line in output.splitlines() if ' = ' in line)
        assert math.isclose(float(metadata['sigma_total_barn']),interpolate(total_grid,total,energy),rel_tol=1e-13)
        assert math.isclose(float(metadata['sigma_elastic_barn']),interpolate(elastic_grid,elastic,energy),rel_tol=1e-13)
        coeff = coefficients(nodes,energy)
        assert math.isclose(float(metadata['mean_mu_cm']),coeff[0],rel_tol=1e-13)
        rows = list(csv.DictReader(io.StringIO('\n'.join(l for l in output.splitlines() if not l.startswith('#')))))
        assert len(rows) == n
        mus = [float(row['mu_cm']) for row in rows]
        mean = sum(mus)/n
        variance = (1+2*coeff[1])/3-coeff[0]**2
        assert abs(mean-coeff[0]) < 6*math.sqrt(variance/n)
        for j in range(10):
            left,right = -1+j*0.2,-1+(j+1)*0.2
            probability = primitive(coeff,right)-primitive(coeff,left)
            count = sum(left <= mu < right for mu in mus)
            assert abs(count-n*probability) < 6*math.sqrt(n*probability*(1-probability))+2
        a = manifest['AWR']
        ratios = [float(row['E_out_eV'])/energy for row in rows]
        expected = (a*a+1+2*a*coeff[0])/(a+1)**2
        sigma_ratio = 2*a/(a+1)**2*math.sqrt(variance/n)
        assert abs(sum(ratios)/n-expected) < 6*sigma_ratio
        quadrants = [0]*4
        for row,mu,ratio in zip(rows,mus,ratios):
            x,y,z = (float(row['omega_'+axis]) for axis in 'xyz')
            assert abs(x*x+y*y+z*z-1) < 1e-12
            assert abs(z-(1+a*mu)/math.sqrt(a*a+1+2*a*mu)) < 1e-12
            assert abs(ratio-(a*a+1+2*a*mu)/(a+1)**2) < 1e-12
            quadrant = int((math.atan2(y,x)%(2*math.pi))/(math.pi/2))
            quadrants[quadrant] += 1
        assert all(abs(count-n/4) < 6*math.sqrt(n*3/16) for count in quadrants)
        print(f'PASS: independent ENDF reference and {n} random collisions at {energy:g} eV; '
              f'<mu_cm>={mean:.6f}, reference={coeff[0]:.6f}')


if __name__ == '__main__':
    main()
