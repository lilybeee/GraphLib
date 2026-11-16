/*
 * png_writer.h — minimal PNG encoder (internal).  8-bit RGB, no interlacing,
 * Sub filter + LZ77/fixed-Huffman deflate, so flat backgrounds compress well.
 */
#ifndef PNG_WRITER_H
#define PNG_WRITER_H

#include <stdio.h>

#include "graph_status.h"

/** Writes @rgb (w*h*3 bytes, row-major, top row first) as a PNG.
 *  GRAPH_ERR_IO on a write error, GRAPH_ERR_OUT_OF_MEMORY, GRAPH_ERR_TOO_LARGE. */
GraphStatus png_write_rgb(FILE *out, const unsigned char *rgb, unsigned w, unsigned h);

#endif /* PNG_WRITER_H */
