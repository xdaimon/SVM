#include "pch.h"
#include "InputHandler.h"
#include "Logger.h"

SVM::InputHandler::InputHandler(ClassifierAdapter* ca) {
	classifierAdapter = ca;
}

void SVM::InputHandler::OnPointerDown(float x, float y, double timeWhenDown) {
	x /= width;
	y /= height;
	pointerDown = true;
	this->timeWhenDown = timeWhenDown;
	PointerCoordinates[0] = x;
	PointerCoordinates[1] = y;
	PointerDownCoord[0] = x;
	PointerDownCoord[1] = y;
}

void SVM::InputHandler::OnPointerMove(float x, float y, double timeWhenMove) {
	if (!pointerDown)
		return;
	x /= width;
	y /= height;
	PointerDx[0] = x - PointerCoordinates[0];
	PointerDx[1] = y - PointerCoordinates[1];
	PointerCoordinates[0] = x;
	PointerCoordinates[1] = y;
	Position2D[0] += PointerDx[0];
	Position2D[1] += PointerDx[1];
	OnUserDrag();
}

void SVM::InputHandler::OnPointerUp(float x, float y, double timeWhenUp, int currentClassLabel, bool is3D) {
	pointerDown = false;
	if (is3D)
		return;
	x /= width;
	y /= height;
	double time = timeWhenUp - timeWhenDown;
	if (time < 1.1) {
		Vector<double> P = Vector<double>(2);
		P.SetComponent(0, x);
		P.SetComponent(1, y);
		Vector<double> Q = Vector<double>(2, PointerDownCoord, false);
		double length = Vector<double>::Length(Vector<double>::Subtract(P, Q));
		if (.0001f < length)
			return;
		double* coords = (double*)malloc(2 * sizeof(double));
		coords[0] = (PointerDownCoord[0] - Position2D[0]) * 2.f - 1.f;
		coords[1] = (PointerDownCoord[1] - Position2D[1]) * 2.f - 1.f;
		coords[0] *= width / height;
		coords[1] *= -1.f;
		coords[0] *= 2.;
		coords[1] *= 2.;
		classifierAdapter->AddTrainExample(Vector<double>(2, coords, true), currentClassLabel);
	}
}

void SVM::InputHandler::SetPointerCoordinateScale(int width, int height) {
	this->width = (float)width;
	this->height = (float)height;
}

void SVM::InputHandler::OnUserDrag() {
	classifierAdapter->OnInputChange(Vector<double>(2, Position2D, false));
}

void SVM::InputHandler::Clear() {
	this->PointerCoordinates[0] = 0.0;
	this->PointerCoordinates[1] = 0.0;
	this->PointerDownCoord[0] = 0.0;
	this->PointerDownCoord[1] = 0.0;
	this->PointerDx[0] = 0.0;
	this->PointerDx[1] = 0.0;
	this->Position2D[0] = 0.0;
	this->Position2D[1] = 0.0;
}
