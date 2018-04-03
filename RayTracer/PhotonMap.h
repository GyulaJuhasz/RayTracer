#ifndef PHOTON_MAP_H
#define PHOTON_MAP_H

#include "Color.h"
#include "Config.h"

namespace raytracer {
	class PhotonMap {
		raytracer::Color cells[FOTONTERKEP_X][FOTONTERKEP_Y];
		bool postProcessDone;

	public:

		PhotonMap();

	private:

		void PostProcess();

	public:

		raytracer::Color GetPhoton(int x, int y);

		void AddPhoton(int x, int y, raytracer::Color photonIntensity);

		void ReadFromFile(const char *fileName);

		void WriteToFile(const char *fileName);
	};
}

#endif // !PHOTON_MAP_H