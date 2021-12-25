#include "PhotonMap.h"
#include <fstream>

using raytracer::PhotonMap;

PhotonMap::PhotonMap() { postProcessDone = false; }

void PhotonMap::PostProcess() {
	for (int a = 0; a < MAX_ATLAG; a++) {
		raytracer::Color temp;
		for (int i = ATLAG_SUGAR; i < FOTONTERKEP_X - ATLAG_SUGAR; i++) {
			for (int j = ATLAG_SUGAR; j < FOTONTERKEP_Y - ATLAG_SUGAR; j++) {
				temp.Set(0.0, 0.0, 0.0);
				for (int x = i - ATLAG_SUGAR; x <= i + ATLAG_SUGAR; x++) {
					for (int y = j - ATLAG_SUGAR; y <= j + ATLAG_SUGAR; y++) {
						raytracer::Color original = cells[x][y];
						temp = temp + cells[x][y];
					}
				}

				temp = temp / ((2 * ATLAG_SUGAR + 1) * (2 * ATLAG_SUGAR + 1));
				cells[i][j] = temp;
			}
		}
	}
	postProcessDone = true;
}

raytracer::Color PhotonMap::GetPhoton(int x, int y) {
	if (!postProcessDone) PostProcess();
	return cells[x][y];
}

void PhotonMap::AddPhoton(int x, int y, raytracer::Color photonIntensity) {
	if ((x >= 0) && (x <= (FOTONTERKEP_X - 1)) && (y >= 0) && (y <= (FOTONTERKEP_Y - 1))) {
		cells[x][y] = cells[x][y] + photonIntensity;
		postProcessDone = false;
	}
}


void PhotonMap::ReadFromFile(const char *fileName) {
	std::ifstream file(fileName);

	if (!file.is_open()) return;

	double R, G, B;

	for (int i = 0; i < FOTONTERKEP_X; i++) {
		for (int j = 0; j < FOTONTERKEP_Y; j++) {
			file >> R >> G >> B;
			cells[i][j].R(R).G(G).B(B);
		}
	}

	file >> postProcessDone;

	file.close();
}

void PhotonMap::WriteToFile(const char *fileName) {
	std::ofstream file(fileName);

	if (!file.is_open()) return;

	for (int i = 0; i < FOTONTERKEP_X; i++) {
		for (int j = 0; j < FOTONTERKEP_Y; j++) {
			file << cells[i][j].GetR() << cells[i][j].GetG() << cells[i][j].GetB();
		}
	}

	file << postProcessDone;

	file.close();
}