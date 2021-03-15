SOURCES = \
    make.jl \
    src/basics.md \
    src/index.md \
    src/install.md \
    src/interpolation.md \
    src/kernels.md \
    src/library.md

default: build/index.html

clean:
	rm -f *~ */*~ */*/*~

build/index.html: $(SOURCES)
	julia --color=yes make.jl

.PHOONY: default clean
