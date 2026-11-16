CC      ?= cc
CFLAGS  ?= -std=c11 -Wall -Wextra -Werror -pedantic
CPPFLAGS += -Iinclude -Isrc -Itest
LDLIBS  += -lm

BUILD := build
OBJDIR := $(BUILD)/obj

HEADERS := $(wildcard include/*.h src/*.h test/*.h)

# The library: everything the unified API in graph.h is built from.
LIB_SRC := $(wildcard src/*.c)
LIB_OBJ := $(LIB_SRC:src/%.c=$(OBJDIR)/%.o)
LIB     := $(BUILD)/libgraph.a

LIBTESTS := analysis connectivity paths flow combinatorial network ops_io util errors advanced
TESTS := list_test matrix_test matrix_backend_test adaptive_test graph_test facade_test graph_advanced_test memory_budget_test graph_social_test graph_stream_test graph_shard_test viz_test showcase \
         $(LIBTESTS:%=list_%_test) $(LIBTESTS:%=matrix_%_test)

.PHONY: all lib test bench crossover demo showcase social examples images gallery million bench-million verify verify-quick clean

all: lib $(TESTS:%=$(BUILD)/%) $(BUILD)/crossover_report $(BUILD)/viz_demo $(BUILD)/social_network

lib: $(LIB)

$(BUILD) $(OBJDIR):
	mkdir -p $@

$(OBJDIR)/%.o: src/%.c $(HEADERS) | $(OBJDIR)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

$(LIB): $(LIB_OBJ) | $(BUILD)
	ar rcs $@ $^

# The two single-representation suites test one module on its own.
$(BUILD)/list_test: test/test_adjacency_list.c $(OBJDIR)/adjacency_list.o $(OBJDIR)/graph_status.o $(OBJDIR)/memory_budget.o $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(OBJDIR)/adjacency_list.o $(OBJDIR)/graph_status.o $(OBJDIR)/memory_budget.o $(LDLIBS) -o $@

$(BUILD)/matrix_test: test/test_adjacency_matrix.c $(OBJDIR)/adjacency_matrix.o $(OBJDIR)/graph_status.o $(OBJDIR)/memory_budget.o $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(OBJDIR)/adjacency_matrix.o $(OBJDIR)/graph_status.o $(OBJDIR)/memory_budget.o $(LDLIBS) -o $@

$(BUILD)/matrix_backend_test: test/test_matrix_backend.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

# adaptive_test uses test/stub_backend.c for exact cost-ratio control.
$(BUILD)/adaptive_test: test/test_adaptive.c test/stub_backend.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) test/test_adaptive.c test/stub_backend.c $(LIB) $(LDLIBS) -o $@

# list_<name>_test / matrix_<name>_test are test/test_list_<name>.c / test_matrix_<name>.c.
$(BUILD)/list_%_test: test/test_list_%.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/matrix_%_test: test/test_matrix_%.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/bench_list: test/bench_list.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/bench_matrix: test/bench_matrix.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

# Timing sweeps (CSV on stdout).  Use CFLAGS="-std=c11 -O2" for real numbers.
bench: $(BUILD)/bench_list $(BUILD)/bench_matrix
	./$(BUILD)/bench_list
	./$(BUILD)/bench_matrix

$(BUILD)/facade_test: test/test_facade.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -pthread -o $@

$(BUILD)/graph_advanced_test: test/test_graph_advanced.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/scale_test: test/test_scale.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/memory_budget_test: test/test_memory_budget.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/graph_social_test: test/test_graph_social.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/graph_stream_test: test/test_graph_stream.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/graph_shard_test: test/test_graph_shard.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/verify_report: test/verify_report.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

# The terminal verification report: expected next to actual, for every area of the library.
# `make verify` runs everything (about a minute); `make verify-quick` skips the million-vertex graphs.
verify: $(BUILD)/verify_report
	./$(BUILD)/verify_report

verify-quick: $(BUILD)/verify_report
	./$(BUILD)/verify_report quick

$(BUILD)/viz_test: test/test_viz.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

$(BUILD)/showcase: test/showcase.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

# Expected vs actual for every algorithm on fixed graphs (writes build/showcase_U*.svg/png).
showcase: $(BUILD)/showcase
	./$(BUILD)/showcase

# The README's example programs, built and run (their output is quoted in the README).
EXAMPLES := quickstart algorithms big_graph advanced sharded readme_images gallery million_users

$(BUILD)/examples/%: examples/%.c $(LIB) $(HEADERS) | $(BUILD)
	mkdir -p $(BUILD)/examples
	$(CC) $(CPPFLAGS) $(CFLAGS) -DBUILD_CFLAGS='"$(CFLAGS)"' $< $(LIB) $(LDLIBS) -o $@

examples: $(EXAMPLES:%=$(BUILD)/examples/%)
	@for e in quickstart algorithms big_graph advanced sharded; do echo "=== $$e"; ./$(BUILD)/examples/$$e; done

# The pictures in docs/images/ (drawn by the library itself; see examples/readme_images.c).
images: $(BUILD)/examples/readme_images $(BUILD)/examples/gallery
	./$(BUILD)/examples/readme_images docs/images
	mkdir -p docs/images/gallery
	./$(BUILD)/examples/gallery docs/images/gallery

# THE BENCHMARK: a million users analysed and drawn, with a report (docs/images/million/REPORT.md and
# results.csv).  The default flags above have no optimisation, so timings from `make` alone are not
# representative; this target rebuilds everything with -O2 in its own directory (build-o2/) first.
# Compare two runs with:  make bench-million BASELINE=old/results.csv
BENCH_CFLAGS := -std=c11 -O2 -Wall -Wextra -Werror -pedantic
bench-million:
	$(MAKE) BUILD=build-o2 CFLAGS="$(BENCH_CFLAGS)" build-o2/examples/million_users
	mkdir -p docs/images/million
	./build-o2/examples/million_users docs/images/million 1000000 $(BASELINE)

# The same example with the default (unoptimised) flags: fast to build, slow to run, not a benchmark.
million: $(BUILD)/examples/million_users
	mkdir -p docs/images/million
	./$(BUILD)/examples/million_users docs/images/million

# One picture per kind and size of graph (docs/images/gallery/); slow part is the 100 000-user network.
gallery: $(BUILD)/examples/gallery
	mkdir -p docs/images/gallery
	./$(BUILD)/examples/gallery docs/images/gallery

$(BUILD)/social_network: examples/social_network.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -pthread -o $@

# A society grows from 6 people to 20000; analysis, pictures in build/social/.  Use -O2:
#   make social CFLAGS="-std=c11 -O2"
social: $(BUILD)/social_network
	mkdir -p $(BUILD)/social
	./$(BUILD)/social_network

$(BUILD)/crossover_report: tools/crossover_report.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

# Prints where the matrix starts beating the list, per algorithm.  Use -O2:
#   make crossover CFLAGS="-std=c11 -O2"
crossover: $(BUILD)/crossover_report
	./$(BUILD)/crossover_report

$(BUILD)/viz_demo: tools/viz_demo.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

# Example pictures in build/demo/ (open the .svg in a browser, the .png in any viewer).
demo: $(BUILD)/viz_demo
	mkdir -p $(BUILD)/demo
	./$(BUILD)/viz_demo random 60 $(BUILD)/demo/random60.svg
	./$(BUILD)/viz_demo clusters 400 $(BUILD)/demo/clusters400.png
	./$(BUILD)/viz_demo scalefree 5000 $(BUILD)/demo/scalefree5000.png
	./$(BUILD)/viz_demo grid 10000 $(BUILD)/demo/grid10000.png
	./$(BUILD)/viz_demo scalefree 100000 $(BUILD)/demo/scalefree100000.png
	./$(BUILD)/viz_demo clusters 2000 $(BUILD)/demo/clusters2000_matrix.png matrix

$(BUILD)/graph_test: test/test_graph.c $(LIB) $(HEADERS) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) $< $(LIB) $(LDLIBS) -o $@

test: all
	@for t in $(TESTS); do echo "== $$t"; ./$(BUILD)/$$t || exit 1; done

clean:
	rm -rf $(BUILD)
