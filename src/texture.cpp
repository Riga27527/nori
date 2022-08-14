
#include <nori/texture.h>

// #define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>
#include <iostream>
NORI_NAMESPACE_BEGIN


// Texture::Texture(const Texture& src) : image_data(src.image_data.clone()), iw(src.iw), 
// 										ih(src.ih), ic(src.ic){ }

Color3f Texture::sampleBilinear(float x, float y) const{
	int fx = (int)x;
	int fy = (int)y;
	int x1 = clamp(fx, 0, iw-1);
	int y1 = clamp(fy, 0, ih-1);
	int x2 = clamp(x1 + 1, 0, iw-1);
	int y2 = clamp(y1 + 1, 0, ih-1);
	float dx = x - (float)fx;
	float dy = y - (float)fy;
	Color3f c00 = getPixel(x1, y1), c01 = getPixel(x2, y1), c10 = getPixel(x1, y2), c11 = getPixel(x2, y2);
	Color3f c1 = c00 * (1.0f - dx) + c01 * dx;
	Color3f c2 = c10 * (1.0f - dx) + c11 * dx;
	Color3f c = c1 * (1.0f - dy) + c2 * dy;
	return c;
}

NORI_NAMESPACE_END