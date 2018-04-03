#ifndef FILE_UTILS_H
#define FILE_UTILS_H
#endif

#include <algorithm>
#include <fstream>

namespace raytracer {

	namespace files {

		inline void saveToBmp(const char* fileName, const float* image, const int imageWidth, const int imageHeight) {
			std::ofstream file(fileName);

			if (file.is_open()) {

				int BMPheadersize = 14;
				int DIBheadersize = 40;

				// One pixel takes 3 bytes
				int RowSize = imageWidth * 3;

				int RowPadding = (4 - ((imageWidth * 3) % 4)) % 4;

				// Total image size in bytes
				int BitmapSize = imageHeight * (RowSize + RowPadding);

				int FileSize = BMPheadersize + DIBheadersize + BitmapSize;

				/// BMP Header start

				// BMP format characters
				file << (unsigned char)'B';
				file << (unsigned char)'M';

				// Size of BMP in little-endian
				for (int i = 0; i < 4; i++) {
					file << (unsigned char)(FileSize % 256);
					FileSize /= 256;
				}

				// 0 for reserved bits
				for (int i = 0; i < 4; i++) {
					file << (unsigned char)0;
				}

				int BitmapOffset = BMPheadersize + DIBheadersize;

				// Image offset in little-endian
				for (int i = 0; i < 4; i++) {
					file << (unsigned char)(BitmapOffset % 256);
					BitmapOffset /= 256;
				}

				// DIB header size in little-endian
				for (int i = 0; i < 4; i++) {
					file << (unsigned char)(DIBheadersize % 256);
					DIBheadersize /= 256;
				}

				int Width = imageWidth;
				int Height = imageHeight;

				// Image width in little-endian
				for (int i = 0; i < 4; i++) {
					file << (unsigned char)(Width % 256);
					Width /= 256;
				}

				// Image height in little-endian
				for (int i = 0; i < 4; i++) {
					file << (unsigned char)(Height % 256);
					Height /= 256;
				}

				// Number of levels = 1  (  01 00 =>  00 01)
				file << (unsigned char)1;
				file << (unsigned char)0;

				// Bits / pixel
				file << (unsigned char)24;
				file << (unsigned char)0;

				// Plain RGB format
				for (int i = 0; i < 4; i++) {
					file << (unsigned char)0;
				}

				// Size of image data
				for (int i = 0; i < 4; i++) {
					file << (unsigned char)(BitmapSize % 256);
					BitmapSize /= 256;
				}

				// Horizontal resolution
				file << (unsigned char)13;
				file << (unsigned char)11;
				file << (unsigned char)0;
				file << (unsigned char)0;

				// Vertical resolution
				file << (unsigned char)13;
				file << (unsigned char)11;
				file << (unsigned char)0;
				file << (unsigned char)0;

				// Number of palette colors
				for (int i = 0; i < 4; i++) {
					file << (unsigned char)0;
				}

				// Number of important palette colors
				for (int i = 0; i < 4; i++) {
					file << (unsigned char)0;
				}

				/// Raw data

				int index = 0;

				double tempR;
				double tempG;
				double tempB;

				int R, G, B;

				for (int i = 0; i < imageHeight; i++) {
					for (int j = 0; j < imageWidth; j++) {

						index = i * 3 * imageWidth + j * 3;

						tempR = image[index] * 255;
						tempG = image[index + 1] * 255;
						tempB = image[index + 2] * 255;

						R = min((int)floor(tempR + 0.5f), 255);
						G = min((int)floor(tempG + 0.5f), 255);
						B = min((int)floor(tempB + 0.5f), 255);

						// B, G, R in order
						file << (unsigned char)B;
						file << (unsigned char)G;
						file << (unsigned char)R;

					}

					// Zeroes for padding
					for (int i = 0; i < RowPadding; i++) {
						file << (unsigned char)0;
					}
				}

				file.close();
			}
		}
	}

}