#pragma once

#include "ClassifierAdapter.h"

namespace SVM
{
	class MSVMAdapter : public ClassifierAdapter
	{
	public:
		MSVMAdapter() { classifier = new MulticlassSVM(); };
		//~MSVMAdapter() { delete classifier; };

		void Setup(SVM::RenderState state);
		void Update();
		void Train();
		void TrainStep(double step_size, int step_itters, double regularization);
		void PrintInformation();
		void Clear();

		void AddTrainExample(Vector<double> const& example, int const label);
		void AddTestExample(Vector<double> const& example, int const label);
		void OnInputChange(Vector<double> positionInfo);

		void SetDataSet(int i);

		RenderingResources& GetRenderingResources();
		Vector<double> GetClassColor(int i);

	private:
		MulticlassSVM* classifier;

		int steps_taken = 0;

		int currentDataSet = 0;

		static const float Colors[];

		void GeneratePlottingResources(RenderingResources &rr, int dimension);
		void GenerateWeightsDisplayResources(RenderingResources &rr);
		std::string Generate3dShader();
		std::string Generate2dShader();
		std::string GenerateWeightsDisplayShader();
		void SetupRenderingResources();
		RenderingResources renderingResources;

		// All following Matrices are stored column major in memory
		Matrix<double> TrainExamples;
		Vector<int> TrainLabels;
		Matrix<double> TestExamples;
		Vector<int> TestLabels;

		// A validation set?
		// Matrix validationExamples;
		// Vector validationLabels;
		void GenerateRandomData(SVM::RenderState state);
		void MakeImageDataPoints();

	};
}
