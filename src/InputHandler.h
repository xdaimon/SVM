#pragma once

#include "ClassifierAdapter.h"
#include "Vector.h"
#include "Matrix.h"

namespace SVM
{
	class InputHandler
	{
	public:
		InputHandler(ClassifierAdapter* ca);

		void OnPointerDown(float x, float y, double timeWhenDown);
		void OnPointerMove(float x, float y, double timeWhenMove);
		void OnPointerUp(float x, float y, double timeWhenUp, int currentClassLabel, bool is3D);

		void OnKeyDown();
		void OnKeyUp();

		void Clear();

		void SetPointerCoordinateScale(int width, int height);

	private:
		void OnUserDrag();

		double PointerCoordinates[2];
		double PointerDownCoord[2];
		double PointerDx[2];
		double Position2D[2];

		ClassifierAdapter* classifierAdapter;

		bool pointerDown;
		unsigned long long timeWhenDown;
		float width;
		float height;
	};
}
