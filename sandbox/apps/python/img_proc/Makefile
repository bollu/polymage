APPS = bilateral_grid \
	   campipe \
	   harris \
	   interpolate \
	   pyramid_blend \
	   unsharp_mask

all: 
	@-for d in $(APPS); do \
		make -C $$d $@; \
		done

clean:
	@-for d in $(APPS); do \
		make -C $$d clean; \
		done
