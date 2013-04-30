# Copyright (c) 2012, Richard P. Curnow
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the <organization> nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

CC := gcc
#CFLAGS := -g -Wall
CFLAGS := -Os -Wall

PROGS := height4 wm_to_grid grid_to_wm create_picture create_picture_2 etrs89
all: $(PROGS)

height4: height4.o coords.o solve.o uk.o contour.o tikz.o svg.o osxx02.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

wm_to_grid: wm_to_grid.o coords.o solve2.o estrin2.o uk.o contour.o tikz.o svg.o osxx02.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

grid_to_wm: grid_to_wm.o coords.o solve2.o estrin2.o uk.o contour.o tikz.o svg.o osxx02.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

create_picture: create_picture.c tikz.o contour.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

create_picture_2: create_picture_2.c tikz.o contour.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

etrs89: etrs89.o coords.o uk.o solve2.o contour.o tikz.o svg.o osxx02.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

distance2: distance2.o solve.o coords.o vincenty.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

distance: distance.o solve.o
	$(CC) $(CFLAGS) -o $@ $^ -lm

%.o : %.c tool.h contour.h
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	-rm -f $(PROGS) *.o

GRID_PNGS := east_new.png north_new.png total_new.png total_10_new.png
pictures: $(GRID_PNGS) height4_new.png

%.png : %.svg
	inkscape --export-png=$@ -D $^

wm_to_grid.out : wm_to_grid
	./wm_to_grid > $@

wm_to_grid_svg.out : wm_to_grid
	./wm_to_grid svg > $@

grid_to_wm.out : grid_to_wm
	./grid_to_wm > $@

grid_to_wm_svg.out : grid_to_wm
	./grid_to_wm svg > $@

height4_new.svg height4.out : height4
	./height4 > height4.out

picture1.tex : create_picture
	./create_picture > $@

picture2.tex : create_picture_2
	./create_picture_2 > $@

pdf: wmen_paper.pdf

wmen_paper.pdf : wmen_paper.tex figures/*.pdf
	pdflatex $<
	pdflatex $<

figs:
	$(MAKE) -j4 -B -f wmen_paper.makefile


