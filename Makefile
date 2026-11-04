FC = gfortran
FFLAGS = -std=f2008 -O2 -Wall -Wextra -Wimplicit-interface -fcheck=all -fbacktrace
BUILD = build

.PHONY: all test clean
all: $(BUILD)/scatter_demo

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/neutron_collision.o: src/neutron_collision.f90 | $(BUILD)
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) -c $< -o $@

$(BUILD)/scatter_demo: app/scatter_demo.f90 $(BUILD)/neutron_collision.o
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) $^ -o $@

$(BUILD)/test_collision: tests/test_collision.f90 $(BUILD)/neutron_collision.o
	$(FC) $(FFLAGS) -J$(BUILD) -I$(BUILD) $^ -o $@

test: $(BUILD)/test_collision $(BUILD)/scatter_demo
	$(BUILD)/test_collision
	python3 tests/check_reference.py
	python3 tests/check_errors.py

clean:
	rm -rf $(BUILD)
