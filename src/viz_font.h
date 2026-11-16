/* viz_font.h — bitmap font used by the PNG renderer (internal). */
#ifndef VIZ_FONT_INCLUDED_H
#define VIZ_FONT_INCLUDED_H

#define VIZ_FONT_W 8
#define VIZ_FONT_H 16

/** VIZ_FONT[c - 32] holds the glyph for ASCII c in 32..126: VIZ_FONT_H rows, MSB = leftmost pixel. */
extern const unsigned char VIZ_FONT[95][VIZ_FONT_H];

#endif /* VIZ_FONT_INCLUDED_H */
