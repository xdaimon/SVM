#include "pch.h"
#include "MulticlassSVM.h"
#include "Matrix.h"
#include "Vector.h"
#include "RenderingResources.h"
#include "Logger.h"
#include <ctime>

using namespace SVM;

void MulticlassSVM::Clear() {
	trainStepItter = 0;
	previous_loss = 1.;
	Weights = Matrix<double>(0, 0, false);
}

void MulticlassSVM::FindLabelsInUse(Vector<int> const& TrainLabels) {
	int label = TrainLabels.GetComponent(0);
	RowToClassMap = Vector<int>(0);
	RowToClassMap.Append(label);
	for (int i = 1; i < TrainLabels.M; ++i) {
		label = TrainLabels.GetComponent(i);
		bool contains_label = false;
		for (int j = 0; j < RowToClassMap.M; ++j)
			contains_label |= label == RowToClassMap.GetComponent(j);
		if (!contains_label)
			RowToClassMap.Append(label);
	}
	for (int i = 0; i < RowToClassMap.M; ++i)
		RowToClassMap.SetComponent(i, RowToClassMap.GetComponent(i) - 1);
}

void MulticlassSVM::TrainStep(Matrix<double> const& TrainExamples, Vector<int> const& TrainLabels, double step_size, int step_itters, double regularization) {
	if (previous_loss < 0.000001)
		return;
	FindLabelsInUse(TrainLabels);
	Vector<int> ClassToRowMap = Vector<int>(SVM::MaxNumClasses);
	for (int i = 0; i < RowToClassMap.M; ++i)
		ClassToRowMap.SetComponent(RowToClassMap.GetComponent(i), i);

	// Initialize Weights
	const int numberOfClasses = RowToClassMap.M;
	std::srand((unsigned)(std::time(0)));
	if (Weights.M * Weights.N == 0 || Weights.M != numberOfClasses) {
		Vector<double> weight_row = Vector<double>(TrainExamples.M);
		Weights = Matrix<double>(numberOfClasses, TrainExamples.M, false);
		for (int i = 0; i < numberOfClasses; ++i) {
			for (int j = 0; j < weight_row.M; ++j)
				// weight_row.SetComponent(j, 0.0 );
				weight_row.SetComponent(j, (std::rand() % 1000) / 10000.0);
			Weights.SetRow(i, weight_row);
		}
	}

	// Return while having the random matrix once
	if (trainStepItter == 0) {
		trainStepItter++;
		return;
	}

	weights_regularization = regularization;
	const double eta = step_size;
	// Optimize loss
	for (int i = 0; i < step_itters; ++i)
		Weights.Add(LossGradient(TrainExamples, TrainLabels, ClassToRowMap).Multiply(-eta));

	// double loss = Loss(TrainExamples, TrainLabels, ClassToRowMap) / (double) TrainExamples.N;
	// LOG("stepSize : %f \tstepItters : %d \tloss : %12.12f", step_size, step_itters, loss);
	// LOG("stepSize : %f \tstepItters : %d ", step_size, step_itters);
}

// random stuff
void MulticlassSVM::Train(Matrix<double> const& TrainExamples, Vector<int> const& TrainLabels) {
	// TODO
	// this function currently causes un responsivness
	return;

	FindLabelsInUse(TrainLabels);
	Vector<int> ClassToRowMap = Vector<int>(SVM::MaxNumClasses);
	for (int i = 0; i < RowToClassMap.M; ++i)
		ClassToRowMap.SetComponent(RowToClassMap.GetComponent(i), i);

	// TODO should initialize to a random matrix!?
	// Initialize Weights
	Weights = Matrix<double>(RowToClassMap.M, TrainExamples.M, false);
	Vector<double> weight_row = Vector<double>(Weights.N);
	std::srand((unsigned)(std::time(0)));
	for (int i = 0; i < Weights.M; ++i) {
		for (int j = 0; j < weight_row.M; ++j)
			// weight_row.SetComponent(j, std::rand() % 1000 / 1000.0 * 2. - 1. );
			weight_row.SetComponent(j, 0.);
		Weights.SetRow(i, weight_row);
	}

	const double eta = 0.001;
	// Optimize loss
	double loss = 1.0;
	while (loss > 0.0001) {
		for (int i = 0; i < 100; ++i)
			Weights.Add(LossGradient(TrainExamples, TrainLabels, ClassToRowMap).Multiply(-eta));
		loss = Loss(TrainExamples, TrainLabels, ClassToRowMap);
	}
}

// TODO it might be cool to do the kNN

Vector<int> MulticlassSVM::Test(Matrix<double> const& TestExamples) {
	// Score
	const Matrix<double> scores = Matrix<double>::Multiply(Weights, TestExamples);

	// Choose highest score for each testExample and store the label in the return vector
	Vector<int> labels = Vector<int>(TestExamples.N);
	double max_score = 0.0;
	int max_index = 0;
	int label = 0;
	for (int j = 0; j < scores.N; ++j) {
		for (int i = 0; i < scores.M; ++i) {
			if (scores.GetElement(i, j) > max_score) {
				max_score = scores.GetElement(i, j);
				max_index = i;
			}
		}
		label = max_index + 1;
		labels.SetComponent(j, label);
		max_score = 0.0;
		max_index = 0;
	}
	return labels;
}

double MulticlassSVM::Loss(Matrix<double> const& trainExamples, Vector<int> const& trainLabels, Vector<int> const& classToRowMap) {
	double loss = 0.0;

	Matrix<double> scores = Matrix<double>::Multiply(Weights, trainExamples);

	auto correctRow = [&](int i) {
		return classToRowMap.GetComponent(-1 + trainLabels.GetComponent(i));
	};
	Matrix<double> correct_scores = Matrix<double>(scores.M, scores.N, true);
	for (int j = 0; j < scores.N; ++j) {
		const double score = scores.GetElement(correctRow(j), j);
		for (int i = 0; i < scores.M; ++i)
			correct_scores.SetElement(i, j, score);
	}

	Matrix<double> margins = Matrix<double>::Subtract(scores, correct_scores);
	const double delta = 1.0;
	for (int i = 0; i < margins.M; ++i)
		for (int j = 0; j < margins.N; ++j)
			margins.SetElement(i, j, margins.GetElement(i, j) + delta);
	for (int k = 0; k < trainExamples.N; ++k)
		margins.SetElement(correctRow(k), k, 0.0);

	auto max = [](double x, double y) {
		if (x > y)
			return x;
		else
			return y;
	};
	for (int i = 0; i < margins.M; ++i)
		for (int j = 0; j < margins.N; ++j)
			loss += max(0.0, margins.GetElement(i, j));

	for (int i = 0; i < Weights.M; ++i)
		for (int j = 0; j < Weights.N; ++j)
			loss += weights_regularization * pow(Weights.GetElement(i, j), 2);

	return loss / (double)trainExamples.N;
}

// define to use hinge loss (max function), otherwise use sigmoid loss
#define USE_HINGE
Matrix<double> MulticlassSVM::LossGradient(Matrix<double> const& trainExamples, Vector<int> const& trainLabels, Vector<int> const& classToRowMap) {
	Matrix<double> dW = Matrix<double>(Weights.M, Weights.N, false);

	const Matrix<double> scores = Matrix<double>::Multiply(Weights, trainExamples);

	auto correctRow = [&](int i) {
		return classToRowMap.GetComponent(-1 + trainLabels.GetComponent(i));
	};
	Matrix<double> correct_scores = Matrix<double>(scores.M, scores.N, true);
	for (int j = 0; j < scores.N; ++j) {
		const double score = scores.GetElement(correctRow(j), j);
		for (int i = 0; i < scores.M; ++i)
			correct_scores.SetElement(i, j, score);
	}

	Matrix<double> margins = Matrix<double>::Subtract(scores, correct_scores);
	const double delta = 1.;
	margins.AddAll(delta);
	for (int k = 0; k < trainExamples.N; ++k)
		margins.SetElement(correctRow(k), k, 0.0);

#ifndef USE_HINGE
	auto sigmoid = [](double x) { return 1. / (1. + exp(-(x))); };
#endif

	for (int k = 0; k < trainExamples.N; ++k) {
		const int correct_row = correctRow(k);
		const Vector<double> example = trainExamples.GetColumn(k);
		int numberMarginsSatisfied = 0;
#ifndef USE_HINGE
		double sig_sum = 0.;
#endif
		for (int i = 0; i < margins.M; ++i) // for all j that're not the correct label
		{
#ifdef USE_HINGE
			if (margins.GetElement(i, k) > 0.0000001) {
				numberMarginsSatisfied++;
				// margins for all correct rows are 0.0
				//if (i != correct_row)
				dW.SetRow(i, Vector<double>::Sum(dW.GetRow(i), example));
			}
#else
			if (i != correct_row) {
				double sp = sigmoid(margins.GetElement(i, k));
				double dsp = sp * (1. - sp);
				sig_sum += dsp;
				dW.SetRow(i, Vector<double>::Sum(dW.GetRow(i), Vector<double>::Multiply(example, dsp)));
			}
#endif
		}
#ifdef USE_HINGE
		dW.SetRow(correct_row, Vector<double>::Subtract(dW.GetRow(correct_row), Vector<double>::Multiply(example, numberMarginsSatisfied)));
#else
		dW.SetRow(correct_row, Vector<double>::Subtract(dW.GetRow(correct_row), Vector<double>::Multiply(example, sig_sum)));
#endif
	}

	// Regularization
	for (int i = 0; i < dW.M; ++i)
		for (int j = 0; j < dW.N; ++j)
			dW.SetElement(i, j, dW.GetElement(i, j) + weights_regularization * Weights.GetElement(i, j));

	return dW;
}
