CUDA_COMPILE = $(NVCC) $(NVCCFLAGS) -lineinfo -x cu -dc -ccbin $(CXX) --compiler-options="$(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CFLAGS) $(CFLAGS)"

.c.o:
	$(CUDA_COMPILE) -o $@ -c $<

.cu.o:
	$(CUDA_COMPILE) -o $@ -c $<

$(gpu_object): $(objects)
	$(NVCC) -dlink -o $@ $^
