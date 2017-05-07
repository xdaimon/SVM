#include "pch.h"

// #include <png++/png.hpp>

#include "RenderingResources.h"
#include "MutliclassSVMAdapter.h"
#include "Logger.h"

void SVM::MSVMAdapter::AddTrainExample(Vector<double> const& example, int const label) {
	if (TrainLabels.M == SVM::MaxNumExamples)
		return;
	Vector<double> v;
	for (int i = 0; i < example.M; ++i)
		v.Append(example.GetComponent(i));
	v.Append(1.0);
	TrainExamples.AddColumn(TrainExamples.N, v);
	TrainLabels.AddComponent(TrainLabels.M, label);
	Update();
}

void SVM::MSVMAdapter::AddTestExample(Vector<double> const& example, int const label) {
	if (TrainLabels.M == SVM::MaxNumExamples)
		return;
	Vector<double> v;
	for (int i = 0; i < example.M; ++i)
		v.Append(example.GetComponent(i));
	v.Append(1.0);
	TestExamples.AddColumn(TestExamples.N, v);
	TestLabels.AddComponent(TestLabels.M, label);
}

void SVM::MSVMAdapter::Setup(RenderState state) {
	GenerateRandomData(state);
	SetupRenderingResources();
}

void SVM::MSVMAdapter::Update() {
	const int dimension = TrainExamples.M;
	if (dimension <= 4 && dimension > 0) {
		// A MSVMAdapter renderingResources is ordered like this when rendering a graph of planes and points:
		// Examples
		// Labels
		// Weights
		// Number of Examples/Weights
		// Colors
		// Coordinates
		double* data;
		if (TrainExamples.M != 0) {
			data = TrainExamples.GetData();
			for (int i = 0; i < TrainExamples.N * dimension; ++i)
				renderingResources.data[0].data[i] = (float)data[i];
		}

		if (TrainLabels.M != 0) {
			int* dataInt = TrainLabels.GetData();
			for (int i = 0; i < TrainLabels.M; ++i)
				((int*)renderingResources.data[1].data)[i] = dataInt[i];
		}

		if (classifier->Weights.M != 0) {
			data = classifier->Weights.GetData();
			const int num_rows = classifier->Weights.M;
			const int num_cols = classifier->Weights.N;
			for (int i = 0; i < num_rows; ++i)
				for (int j = 0; j < num_cols; ++j)
					renderingResources.data[2].data[i * num_cols + j] = (float)data[i * num_cols + j];
			for (int i = 0; i < num_rows; ++i)
				(renderingResources.data[2].data + SVM::MaxNumClasses * dimension)[i] = (float)classifier->RowToClassMap.GetComponent(i);
		}

		if (renderingResources.data.size() > 1) {
			renderingResources.data[3].data[0] = (float)classifier->Weights.M;
			renderingResources.data[3].data[1] = (float)TrainExamples.N;
		}
	} else {
		double* data;
		if (classifier->Weights.M != 0) {
			data = classifier->Weights.GetData();
			const int num_images = classifier->Weights.M;
			const int num_cols = classifier->Weights.N;
			if (renderingResources.data[0].data != nullptr) {
				for (int img = 0; img < num_images; ++img)
					for (int i = 0; i < 180; ++i)
						for (int j = 0; j < 300; ++j)
							renderingResources.data[0].data[img * 256 * 512 + i * 512 + j] = (float)data[img * num_cols + i * 300 + j];
			}
		}
	}
}

void SVM::MSVMAdapter::Train() {
	if (TrainExamples.M != 0) {
		classifier->Train(TrainExamples, TrainLabels);
	}
	Update();
}

void SVM::MSVMAdapter::TrainStep(double step_size, int step_itters, double regularization) {
	if (TrainExamples.M != 0) {
		classifier->TrainStep(TrainExamples, TrainLabels, step_size, step_itters, regularization);
		steps_taken += step_itters;
	}
	Update();
}

void SVM::MSVMAdapter::PrintInformation() {
	// LOG("Number of data points: %d", TrainExamples.N);
	// double accuracy = 0.0;
	// Vector<int> predicted_labels = classifier->Test(TrainExamples);
	// for (int i = 0; i < predicted_labels.M; ++i)
	// 	if (predicted_labels.GetComponent(i) == TrainLabels.GetComponent(i))
	// 		accuracy+=1.0;
	// accuracy /= TrainLabels.M;
	// LOG("Test accuracy on train set: %f", accuracy);
	LOG("Weights:");
	for (int i = 0; i < classifier->Weights.M; ++i) {
		for (int j = 0; j < classifier->Weights.N; ++j)
			printf("%f\t", classifier->Weights.GetElement(i, j));
		printf("\n");
		printf("\n");
	}
	LOG("Steps taken: %d", steps_taken);
	LOG("");
}

void SVM::MSVMAdapter::Clear() {
	// memset(renderingResources.data[2].data, 0, SVM::MaxNumClasses * TrainExamples.M * sizeof(float));
	TrainExamples = Matrix<double>(0, 0, true);
	TestExamples = Matrix<double>(0, 0, true);
	TrainLabels = Vector<int>(0);
	TestLabels = Vector<int>(0);
	steps_taken = 0;
	classifier->Clear();
	Update();
}

void SVM::MSVMAdapter::OnInputChange(Vector<double> positionInfo) {
	double* data = positionInfo.GetData();
	if (data == nullptr) {
		LOG("Nullptr in OnInputChange");
		exit(-1);
	}
	(renderingResources.data.end() - 1)->data[0] = (float)data[0];
	(renderingResources.data.end() - 1)->data[1] = (float)data[1];
}

using std::string;
using std::to_string;
string SVM::MSVMAdapter::Generate3dShader() {
	string ret = string(R"(
		const int MaxNumClasses = )") +
	             to_string(SVM::MaxNumClasses) + string(R"( ;
		const int MaxNumExamples = )") +
	             to_string(SVM::MaxNumExamples) + string(R"( ;

		uniform vec4 Examples[ MaxNumExamples ];
		uniform int Labels[ MaxNumExamples ];
		uniform vec4 Weights[ MaxNumClasses ];
		uniform vec3 Colors[ MaxNumClasses ];
		uniform float weightToClassMap[ MaxNumClasses ];

		uniform vec2 translation;
		uniform int numberOfClasses;
		uniform int numberOfExamples;

		const float Gamma = 1.81;
		const float Bound = 2.5;
		const float Lens = 2.2;

		const float Shininess = 23.;
		const float PlaneOpacity = .80;
		const float ExampleOpacity = .86;
		const float ExampleSphereRadius = 0.1;

		//#define COLORFUL_ISOLINES

		vec3 Light0;
		vec3 Light1;
		const vec3 SKY = vec3(.85);

		#define EXAMPLE 0
		#define PLANE 1

		out vec4 fragColor;

		struct Intersection {
			// What type of object did we intersect
			int objectID;

			// Which two weights define the plane we've intersected
			ivec2 ids;

			float t;

			float opacity;
		};

		void makeRay(out vec3 ro, out vec3 rd) {
			vec2 pixelCoord = gl_FragCoord.xy / resolution.xy;
			pixelCoord = pixelCoord * 2.0 - 1.0;
			pixelCoord.x *= resolution.x/resolution.y;

			rd = vec3(pixelCoord, Lens);
			rd = normalize(rd);

			ro = vec3(0.0, 0.0, -7.0);

			float cx = cos(10.*translation.x);
			float cy = cos(10.*translation.y);
			float sx = sin(10.*translation.x);
			float sy = sin(-10.*translation.y);
			mat3 rotate = mat3( cx, -sx*sy, sx*cy,
							    0.,     cy,    sy,
							   -sx, -cx*sy, cx*cy );
			rd *= rotate;
			ro *= rotate;
		}

		// Plane-Ray intersection
		float plane(in vec3 ro, in vec3 rd, in vec3 norm, in vec3 pointInPlane) {
			return dot(pointInPlane-ro,norm)/dot(rd, norm);
		}

		// Sphere-Ray intersection
		vec3 sphere(in vec3 ro, in vec3 rd, in vec3 center, in float radius) {
			vec3 oc = ro-center;
			float b = dot( oc, rd );
			float c = dot( oc, oc ) - radius*radius;
			float h = b*b - c;
			return vec3(-b-sqrt(max(h,0.)), -b+sqrt(max(h,0.)), h);
		}

		// TODO inline?
		int maxWeightAtPoint(in vec3 p) {
			// Find the class of the point
			int max_weight = -1;
			float score = 0.;
			float max_score = -1e7;
			for (int i = 0; i < numberOfClasses; ++i) {
				score = dot(vec4(p,1.0), Weights[i]);
				if (score > max_score) {
					max_score = score;
					max_weight = i;
				}
			}
			return max_weight; 
		}

		#define withinBound(t) (t>BoundInterval.x&&t<BoundInterval.y)

		// More than enough intersections for all the planes and examples
		Intersection intersections[MaxNumClasses*(MaxNumClasses-1)/2];
		int numIntersections = 0;

		// Fill the intersections array
		bool intersect(in vec3 ro, in vec3 rd) {
			// Intersect bounding sphere
			vec2 BoundInterval;
			vec3 sphereIntersect = sphere(ro, rd, vec3(0.), Bound);
			// float nearestSphereRayDist = sqrt(max(Bound*Bound-sphereIntersect.z, 0.0))-Bound;
			// if (nearestSphereRayDist < 0.01)
			if (sphereIntersect.z > 0.) {
				// intersections[numIntersections].objectID = BOUND;
				// // intersections[numIntersections].opacity = 1.-clamp(nearestSphereRayDist/0.02, 0., 1.);
				// intersections[numIntersections].opacity = abs(sphereIntersect.y - sphereIntersect.x);
				// intersections[numIntersections].t = sphereIntersect.x;
				// numIntersections++;
				BoundInterval = sphereIntersect.xy;
			} else {
				return false;
			}

			// Find intersections with examples
			for (int i = 0; i < 60; ++i) {
				sphereIntersect = sphere(ro, rd, Examples[i].xyz, ExampleSphereRadius);

				// float nearestSphereRayDist = sqrt(max(ExampleSphereRadius*ExampleSphereRadius-sphereIntersect.z, 0.0))-ExampleSphereRadius;
				// const float closeEnough = ExampleSphereRadius/5.;
				//if (nearestSphereRayDist < closeEnough)

				if (sphereIntersect.z  > 0.) {
					// checking for positive is needed when tracing the shadows??
					if (withinBound(sphereIntersect.x) && sphereIntersect.x > 0.) {
						intersections[numIntersections].objectID = EXAMPLE;

						intersections[numIntersections].ids.x = i;
						intersections[numIntersections].t = sphereIntersect.x;

						// intersections[numIntersections].opacity = 1.-clamp(nearestSphereRayDist/closeEnough, 0., 1.);
						// intersections[numIntersections].opacity *= ExampleOpacity;
						intersections[numIntersections].opacity = clamp(ExampleOpacity * (sphereIntersect.y-sphereIntersect.x)/ExampleSphereRadius, 0., 1.);

						numIntersections++;
					}
				}
			}

			// Find intersections with planes
			//
			// f(x,y,z) is a coefficient matrix for a 3d function formed by the difference between two weight functions
			// The solutions to f(x,y,z) == 0 form a plane. The plane can be parameterized by solving for z and forming
			// the funtion h(x,y). Similarly for two planes the level set of their difference h0(x,y) - h1(x,y) == 0
			// forms a line which can be parameterized by solving for y and forming the function g(x). Given this
			// function g(x) we can perform antialiasing on the seams where two planes intersect
			vec4 f;
			vec3 P;
			float t;
			int max_weight;
			// Find points on ray where any two scores are equivalent
			for (int i = 0; i < numberOfClasses; ++i) {
				for (int j = numberOfClasses-1; j > i; --j) {
					// Construct the coefficients for the linear function f(x,y,z)
					f = Weights[i]-Weights[j];
					// Let P be the z intercept of f(x,y,z) == 0
					P = vec3(0.,0.,-f.w/f.z);
					// Solve for t such that ro+rd*t is on the plane. normailzation of f.xyz (plane normal) cancels out
					t = dot(P-ro, f.xyz)/dot(rd, f.xyz);
					// We're only interested in intersections where the two score functions are maximum.
					// checking for positive is needed when tracing the shadows??
					if (withinBound(t) && t > 0.) {
						// Try to ensure that intersections[k].ids.x holds the id of the max weight in the space before
						// the ray hits the plane -> so find the maxWeight slightly before the intersection
						max_weight = maxWeightAtPoint(ro+rd*(t-0.0001));
						if (max_weight == j || max_weight == i) {
							intersections[numIntersections].objectID = PLANE;
							if (max_weight == j) {
								intersections[numIntersections].ids.x=j;
								intersections[numIntersections].ids.y=i;
							} else {
								intersections[numIntersections].ids.x=i;
								intersections[numIntersections].ids.y=j;
							}
							intersections[numIntersections].t = t;
							intersections[numIntersections].opacity = PlaneOpacity;
							numIntersections++;
						}
					}
				}
			}
			// Sort intersections. Intersection closest to ro is in intersections[0]
			for( int i=0; i<numIntersections-1; i++ )
			for( int j=i+1; j<numIntersections; j++ ) {
				if(intersections[j].t<intersections[i].t) {
					Intersection temp = intersections[i];
					intersections[i] = intersections[j];
					intersections[j] = temp;
				}
			}

			return numIntersections > 0;
		}

		// This shader really needs to take into consideration shadowing. Especially since the specular component
		// shows up in places where there should be no light from the sun.
		vec3 lighting(in vec3 hit, in vec3 rd, in vec3 norm, in vec3 diffuseColor) {
			vec3 ret;

			vec3 lr = normalize(hit-Light0);
			float diff = max(dot(norm, -lr),0.0);
			vec3 refL = reflect(lr, norm);
			float spec = pow(max(dot(refL, -rd),0.0), Shininess);
			vec3 diffComp = diffuseColor * diff;
			vec3 specComp = vec3(spec);
			vec3 ambientComp = diffuseColor*.3;
			
			// Add a flashlight with the camera
			diffComp = mix(diffComp, abs(dot(rd, norm))*diffuseColor, .5);
			ret = diffComp + .2 * specComp + ambientComp;

			lr = normalize(hit - Light1);
			diff = max(dot(norm, -lr), 0.0);
			refL = reflect(lr, norm);
			spec = pow(max(dot(refL, -rd),0.0), Shininess);
			diffComp = diffuseColor * diff;
			specComp = vec3(spec);

			diffComp = mix(diffComp, abs(dot(rd, norm))*diffuseColor, .5);
			ret = mix(diffComp + .2 * specComp + ambientComp, ret, .5);

			return ret;
		}

		vec3 getPlaneColor(in vec3 ro, in vec3 rd, in float t, in int i, in int j) {
			vec3 f = (Weights[i]-Weights[j]).xyz;
			f = normalize(f);
			return lighting(ro + rd*t, rd, f.xyz, Colors[int(weightToClassMap[i])]);
		}

		// k -> kth intersection
		vec3 isolines(in vec3 ro, in vec3 rd, in int k, in vec3 baseColor) {
			// AntiAlias the seams where two planes intersect
			// The goal is to find two points on this seam. A line is constructed from these points.
			// The minimum distance between the primary rays intersecting point Q and
			// the line is found and used for color mixing.
			float dist_to_seam = 1e10;
			vec4 f0;
			vec4 f1;
			vec3 g;
			vec3 z0;
			vec3 z1;
			int idOfPlane;
			float s;
			float t;
			vec3 Q;
			vec2 line;
			vec3 lineP0;
			vec3 lineP1;
			vec3 lineDir;
			Q = ro+rd*intersections[k].t;
			f0 = Weights[intersections[k].ids.x]-Weights[intersections[k].ids.y];
			// A point on the z0 plane is then vec3(s, t, dot(z0, vec3(s,t,1.)))
			z0 = -vec3(f0.x,f0.y,f0.w)/f0.z;
			for (int i = 0; i < numberOfClasses; ++i) {
				//TODO project into screen plane and then take distances in order to smooth out the pixel mixing gradients

				if (i == intersections[k].ids.y)
					continue;

				f1 = Weights[intersections[k].ids.x] - Weights[i];

				// Coefficients for equation of plane
				z1 = -vec3(f1.x,f1.y,f1.w)/f1.z;

				g = z0-z1;

				// A point on the line is vec2(s, dot(line, vec2(s, 1.)))
				line = -vec2(g.x,g.z)/g.y;

				s = 0.;
				t = dot(line, vec2(s, 1.));
				lineP0 = vec3(s, t, dot(z0, vec3(s, t, 1.)));

				s = 1.;
				t = dot(line, vec2(s, 1.));
				lineP1 = vec3(s, t, dot(z0, vec3(s, t, 1.)));

				lineDir = normalize(lineP1-lineP0);

				// Point on line such that point-Q is orthogonal to line
				lineP1 = lineP0+lineDir*dot(lineDir, Q-lineP0);
				float D = length(Q-lineP1);
				if (D < dist_to_seam) {
					dist_to_seam = D;
					idOfPlane = i;
				}
			}

			vec3 colToLerp = getPlaneColor(ro, rd, intersections[k].t, intersections[k].ids.x, idOfPlane);
			#ifdef COLORFUL_ISOLINES
				// return mix(baseColor, colToLerp, .5-.5*smoothstep(0., 20./resolution.x, dist_to_seam));
				return mix(baseColor, colToLerp, .5-.5*smoothstep(0., .022, dist_to_seam));
			#else
				// return mix(baseColor, vec3(0.), .5-.5*smoothstep(0., 20./resolution.x, dist_to_seam));
				return mix(baseColor, .2*colToLerp, .5-.5*smoothstep(0., .022, dist_to_seam));
			#endif
		}

		vec3 getColor(in vec3 ro, in vec3 rd) {
			vec3 color = SKY;
			vec3 objColor = vec3(0.);
			for (int i = numIntersections-1; i >= 0; --i) {
				switch(intersections[i].objectID) {
				case PLANE:
					objColor = getPlaneColor(ro, rd, intersections[i].t, intersections[i].ids.x, intersections[i].ids.y);
					objColor = isolines(ro, rd, i, objColor);
					break;
				case EXAMPLE:
					vec3 hit = ro + rd * intersections[i].t;
					objColor = lighting(hit, rd, normalize(hit-Examples[intersections[i].ids.x].xyz), Colors[Labels[intersections[i].ids.x]-1]);
					break;
				}
				color = mix(color, objColor, intersections[i].opacity);
			}
			return color;
		}
		void main() {
			vec3 ro;
			vec3 rd;
			makeRay(ro, rd);

			Light0 = 1.5*vec3(Bound, .75, Bound);
			Light1 = 1.5*vec3(-Bound, .75, -Bound);

			fragColor.rgb = SKY;
			if ( intersect(ro, rd) ) {
				vec3 color = getColor(ro, rd);
				// if (intersect(toward light))
				// {
				// 	color = shadow;
				// }
				fragColor.rgb = color;
			}
			// TODO intersect another ray from hit in the direction toward light
			// if there is no occluders then intersect() returns false

			// I dont really understand gamma i guess.
			fragColor.rgb = pow(fragColor.rgb, vec3(1./Gamma));
			fragColor.a = 1.;
		}
	)");
	return ret;
}

string SVM::MSVMAdapter::Generate2dShader() {
	string ret = std::string(R"(
		uniform vec3 Examples[)") + to_string(SVM::MaxNumExamples) + string(R"(];
		uniform int Labels[   )") + to_string(SVM::MaxNumExamples) + string(R"(];
		uniform vec3 Weights[ )") + to_string(SVM::MaxNumClasses) + string(R"(];
		uniform vec3 Colors[  )") + to_string(SVM::MaxNumClasses) + string(R"(];
		uniform vec2 translation;
		uniform int numberOfClasses;
		uniform int numberOfExamples;
		uniform float weightToClassMap[)") + to_string(SVM::MaxNumClasses) + string(R"(];

		const float Gamma = 1.5;

		out vec4 fragColor;
		void main() {
			vec2 p = vec2(gl_FragCoord.xy/resolution.xy);
			p -= vec2(translation.x, -translation.y);
			p = p * 2.0 - 1.0;
			p.x *= resolution.x / resolution.y;
			p *= 2.0;

			// Classify position p based on weights
			float score;
			float last_max_score = -1e7;
			int max_i;
			for (int i = 0; i < numberOfClasses; ++i) {
				score = dot(Weights[i], vec3(p, 1.0));
				if (score > last_max_score) {
					last_max_score = score;
					max_i = i;
				}
			}
			float dist_to_seam = 1e10;
			{
				// This is an adaptation of stuff in the 3d shader
				vec3 f;
				float s;
				vec2 line;
				vec2 lineP0;
				vec2 lineP1;
				vec2 lineDir;
				vec2 Q = p;
				for (int i = 0; i < numberOfClasses; ++i) {
					f = Weights[max_i]-Weights[i];
					line = -vec2(f.x, f.z)/f.y;
					s = 0.;
					lineP0 = vec2(s, dot(line, vec2(s, 1.)));
					s = 1.;
					lineP1 = vec2(s, dot(line, vec2(s, 1.)));
					lineDir = normalize(lineP1 - lineP0);
					lineP1 = lineP0+lineDir*dot(lineDir, Q-lineP0);
					float D = length(Q-lineP1);
					if (D < dist_to_seam)
						dist_to_seam = D;
				}
			}
			// Color region
			vec3 color = Colors[int(weightToClassMap[max_i])];
			color *= smoothstep(0.0, 0.02, dist_to_seam);
			// Draw points
			for (int i = 0; i < numberOfExamples; ++i) {
				// color lerp -> 1 when length -> 0
				float color_lerp = 1.0 - smoothstep(0.0, 0.05, length(p - Examples[i].xy));
				float dcolor_lerp = .5*fwidth(color_lerp);
				color = mix(color, vec3(0.0), smoothstep(0.3 - dcolor_lerp,0.45 + dcolor_lerp, color_lerp));
				color = mix(color, Colors[Labels[i] - 1], smoothstep(0.7,0.85, color_lerp));
			}
			fragColor.rgb = color;
			fragColor = pow(fragColor, vec4(1./Gamma));
			fragColor.a = 1.0;
		}
	)");

	return ret;
}

string SVM::MSVMAdapter::GenerateWeightsDisplayShader() {
	string ret = std::string(R"(
		out vec4 fragColor;
		uniform sampler2D weights;
		void main() {
			vec2 uv = gl_FragCoord.xy/resolution.xy;
			uv.x *= resolution.x/resolution.y;
			float tmp = uv.y;
			uv.y = uv.x;
			uv.x = 1.-tmp;

			float intensity = texture(weights, uv).r;

			fragColor = vec4(10.)*clamp(intensity, 0., 1.);
			fragColor.a = 1.;
		}
	)");

	return ret;
}

void SVM::MSVMAdapter::GeneratePlottingResources(RenderingResources& rr, int dimension) {
	ResourceArray ra;

	// Examples
	ra.data = (float*)malloc(SVM::MaxNumExamples * dimension * sizeof(float));
	ra.ownsData = true;
	double* data = TrainExamples.GetData();
	for (int i = 0; i < TrainExamples.N * dimension; ++i)
		ra.data[i] = (float)data[i];
	ra.initializer = [dimension](ResourceArray& ra) {
		GLint program;
		glGetIntegerv(GL_CURRENT_PROGRAM, &program);
		GLint ul = glGetUniformLocation(program, "Examples");
		if (dimension == 4)
			ra.uploader = [=](ResourceArray const& ra) {
				glUniform4fv(ul, SVM::MaxNumExamples, ra.data);
			};
		else
			ra.uploader = [=](ResourceArray const& ra) {
				glUniform3fv(ul, SVM::MaxNumExamples, ra.data);
			};
	};
	rr.data.push_back(std::move(ra));

	// Labels
	ra.data = (float*)malloc(SVM::MaxNumExamples * sizeof(int));
	ra.ownsData = true;
	int* dataInt = TrainLabels.GetData();
	for (int i = 0; i < TrainLabels.M; ++i)
		((int*)ra.data)[i] = dataInt[i];
	ra.initializer = [](ResourceArray& ra) {
		GLint program;
		glGetIntegerv(GL_CURRENT_PROGRAM, &program);
		GLint ul = glGetUniformLocation(program, "Labels");
		ra.uploader = [=](ResourceArray const& ra) {
			glUniform1iv(ul, SVM::MaxNumExamples, (int*)ra.data);
		};
	};
	rr.data.push_back(std::move(ra));

	// Weights and class map
	ra.data = (float*)malloc((SVM::MaxNumClasses * dimension + SVM::MaxNumClasses) * sizeof(float));
	ra.ownsData = true;
	const int num_rows = classifier->Weights.M;
	const int num_cols = classifier->Weights.N;
	if (num_rows * num_cols != 0) {
		data = classifier->Weights.GetData();
		for (int i = 0; i < num_rows; ++i)
			for (int j = 0; j < num_cols; ++j)
				ra.data[i * num_cols + j] = (float)data[i * num_cols + j];
		for (int i = 0; i < num_rows; ++i)
			(ra.data + SVM::MaxNumClasses * dimension)[i] = (float)classifier->RowToClassMap.GetComponent(i);
	}
	ra.initializer = [dimension](ResourceArray& ra) {
		GLint program;
		glGetIntegerv(GL_CURRENT_PROGRAM, &program);
		GLint ul = glGetUniformLocation(program, "Weights");
		GLint ul2 = glGetUniformLocation(program, "weightToClassMap");
		if (dimension == 4)
			ra.uploader = [=](ResourceArray const& ra) {
				glUniform4fv(ul, SVM::MaxNumClasses, ra.data);
				glUniform1fv(ul2, SVM::MaxNumClasses, ra.data + SVM::MaxNumClasses * dimension);
			};
		else
			ra.uploader = [=](ResourceArray const& ra) {
				glUniform3fv(ul, SVM::MaxNumClasses, ra.data);
				glUniform1fv(ul2, SVM::MaxNumClasses, ra.data + SVM::MaxNumClasses * dimension);
			};

	};
	rr.data.push_back(std::move(ra));

	// Number of classes/examples
	ra.data = (float*)malloc(2 * sizeof(float));
	ra.ownsData = true;
	ra.data[0] = (float)classifier->Weights.M;
	ra.data[1] = (float)TrainExamples.N;
	ra.initializer = [](ResourceArray& ra) {
		GLint program;
		glGetIntegerv(GL_CURRENT_PROGRAM, &program);
		GLint ul = glGetUniformLocation(program, "numberOfClasses");
		GLint ul2 = glGetUniformLocation(program, "numberOfExamples");
		ra.uploader = [=](ResourceArray const& ra) {
			glUniform1i(ul, (int)ra.data[0]);
			glUniform1i(ul2, (int)ra.data[1]);
		};
	};
	rr.data.push_back(std::move(ra));

	// Colors
	ra.data = (float*)Colors;
	ra.ownsData = false;
	ra.initializer = [](ResourceArray& ra) {
		GLint program;
		glGetIntegerv(GL_CURRENT_PROGRAM, &program);
		GLint ul = glGetUniformLocation(program, "Colors");
		ra.uploader = [=](ResourceArray const& ra) {
			glUniform3fv(ul, SVM::MaxNumClasses, ra.data);
		};
	};
	rr.data.push_back(std::move(ra));

	// Coordinate data
	ra.data = (float*)malloc(2 * sizeof(float));
	ra.ownsData = true;
	ra.data[0] = ra.data[1] = 0.0f;
	ra.initializer = [](ResourceArray& ra) {
		GLint program;
		glGetIntegerv(GL_CURRENT_PROGRAM, &program);
		GLint ul = glGetUniformLocation(program, "translation");
		ra.uploader = [=](ResourceArray const& ra) {
			glUniform2f(ul, ra.data[0], ra.data[1]);
		};
	};
	rr.data.push_back(std::move(ra));
}

void SVM::MSVMAdapter::GenerateWeightsDisplayResources(RenderingResources& rr) {
	ResourceArray ra;

	// ra.data = (float*)calloc(2 * TrainExamples.M * sizeof(float), 1);
	ra.data = (float*)calloc(512 * 512 * sizeof(float), 1);
	const int num_images = classifier->Weights.M;
	const int num_cols = classifier->Weights.N;
	ra.ownsData = true;
	ra.initializer = [num_images](ResourceArray& ra) {
		GLuint tex;
		glGenTextures(1, &tex);
		glBindTexture(GL_TEXTURE_2D, tex);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, 512, 512, 0, GL_RED, GL_FLOAT, nullptr);
		// glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, num_images*180, 300, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
		ra.uploader = [=](ResourceArray const& ra) {
			glBindTexture(GL_TEXTURE_2D, tex);
			glTexSubImage2D(GL_TEXTURE_2D,
			                0, 0, 0,
			                512,
			                512,
			                GL_RED,
			                GL_FLOAT,
			                ra.data);
		};
	};
	rr.data.push_back(std::move(ra));
}

SVM::RenderingResources& SVM::MSVMAdapter::GetRenderingResources() {
	return renderingResources;
}

Vector<double> SVM::MSVMAdapter::GetClassColor(int index) {
	Vector<double> v = Vector<double>(3);
	v.SetComponent(0, SVM::MSVMAdapter::Colors[index * 3 + 0]);
	v.SetComponent(1, SVM::MSVMAdapter::Colors[index * 3 + 1]);
	v.SetComponent(2, SVM::MSVMAdapter::Colors[index * 3 + 2]);
	return v;
}

// clang-format off
const float SVM::MSVMAdapter::Colors[SVM::MaxNumClasses * 3] = {
	 80 / 255.0f,
	255 / 255.0f,
	  0 / 255.0f,

	  0 / 255.0f,
	131 / 255.0f,
	255 / 255.0f,

	255 / 255.0f,
	220 / 255.0f,
	  0 / 255.0f,

	255 / 255.0f,
	  0 / 255.0f,
	  0 / 255.0f,

	225 / 255.0f,
	225 / 255.0f,
	225 / 255.0f,

	220 / 255.0f,
	  0 / 255.0f,
	120 / 255.0f
}; // clang-format on

void SVM::MSVMAdapter::SetupRenderingResources() {
	RenderingResources rr;
	rr.fragmentShader = string(R"(
		#version 300 es
		precision highp float;
		uniform vec2 resolution;
	 )");

	const int dimension = TrainExamples.M;
	if (dimension == 3 + 1) {
		rr.fragmentShader += Generate3dShader();
		GeneratePlottingResources(rr, 3 + 1);
	} else if (dimension == 2 + 1) {
		rr.fragmentShader += Generate2dShader();
		GeneratePlottingResources(rr, 2 + 1);
	} else {
		rr.fragmentShader += GenerateWeightsDisplayShader();
		GenerateWeightsDisplayResources(rr);
	}
	renderingResources = std::move(rr);
}

// TODO I could generate a set of random points to generate more random points around.
void SVM::MSVMAdapter::GenerateRandomData(RenderState state) {
	const int numberOfDimensions = state == RenderState::Render3D ? 3 : 2;

	std::srand(time(0));
	auto randomReal = [](float lt, float shift) -> double {
		return std::rand() / (double)RAND_MAX * lt + shift;
	};
	auto randomTheta = [randomReal]() -> double {
		return randomReal(3.1415 * 2., 0.0);
	};

	int numExamples;
	if (currentDataSet == 0) {
		numExamples = 60;

		TrainExamples = Matrix<double>(numberOfDimensions + 1, numExamples, true);
		TestExamples = Matrix<double>(numberOfDimensions + 1, numExamples, true);
		TrainLabels = Vector<int>(numExamples);

		for (int i = 0; i < numExamples / 4; ++i) {
			const double th = randomTheta();
			TrainExamples.SetElement(0, i, randomReal(.5, .5) * cos(th));
			TrainExamples.SetElement(1, i, randomReal(.5, .5) * sin(th));
			TrainLabels.SetComponent(i, 1);
		}
		for (int i = numExamples / 4; i < numExamples * 2 / 4; ++i) {
			const double th = randomTheta();
			TrainExamples.SetElement(0, i, randomReal(.25, .0) * cos(th));
			TrainExamples.SetElement(1, i, randomReal(.25, .0) * sin(th));
			TrainLabels.SetComponent(i, 4);
		}
		for (int i = numExamples * 2 / 4; i < numExamples * 3 / 4; ++i) {
			const double th = randomTheta();
			TrainExamples.SetElement(0, i, randomReal(-.25, .3) * cos(th));
			TrainExamples.SetElement(1, i, randomReal(-.75, .1) * sin(th));
			TrainLabels.SetComponent(i, 5);
		}
		for (int i = numExamples * 3 / 4; i < numExamples; ++i) {
			const double th = randomTheta();
			TrainExamples.SetElement(0, i, randomReal(1.05, .7) * cos(th) + .1);
			TrainExamples.SetElement(1, i, randomReal(-.65, .0) * sin(th) - .7);
			TrainLabels.SetComponent(i, 2);
		}
		if (state == RenderState::Render3D) {
			for (int i = 0; i < numExamples / 4; ++i)
				TrainExamples.SetElement(2, i, 1.3);
			for (int i = numExamples / 4; i < numExamples * 2 / 4; ++i)
				TrainExamples.SetElement(2, i, -1.3);
			for (int i = numExamples * 2 / 4; i < numExamples * 3 / 4; ++i)
				TrainExamples.SetElement(2, i, 0.);
			for (int i = numExamples * 3 / 4; i < numExamples; ++i)
				TrainExamples.SetElement(2, i, .42);
		}
	} else {
		numExamples = 60;
		TrainExamples = Matrix<double>(numberOfDimensions + 1, numExamples, true);
		TestExamples = Matrix<double>(numberOfDimensions + 1, numExamples, true);
		TrainLabels = Vector<int>(numExamples);

		for (int i = 0; i < numExamples * 1 / 6; ++i) {
			TrainExamples.SetElement(0, i, randomReal(.4, -.8));
			TrainExamples.SetElement(1, i, randomReal(.4, -.8));
			TrainLabels.SetComponent(i, 1);
		}
		for (int i = numExamples * 1 / 6; i < numExamples * 2 / 6; ++i) {
			TrainExamples.SetElement(0, i, randomReal(.7, .5));
			TrainExamples.SetElement(1, i, randomReal(.7, .5));
			TrainLabels.SetComponent(i, 2);
		}
		for (int i = numExamples * 2 / 6; i < numExamples * 3 / 6; ++i) {
			TrainExamples.SetElement(0, i, randomReal(-0.6, -1.0));
			TrainExamples.SetElement(1, i, randomReal(0.7, 0.5));
			TrainLabels.SetComponent(i, 3);
		}
		for (int i = numExamples * 3 / 6; i < numExamples * 4 / 6; ++i) {
			TrainExamples.SetElement(0, i, randomReal(-0.5, -1.0));
			TrainExamples.SetElement(1, i, randomReal(-0.7, -0.7));
			TrainLabels.SetComponent(i, 4);
		}
		for (int i = numExamples * 4 / 6; i < numExamples * 5 / 6; ++i) {
			TrainExamples.SetElement(0, i, randomReal(1.0, 1.));
			TrainExamples.SetElement(1, i, randomReal(-.8, -.5));
			TrainLabels.SetComponent(i, 5);
		}
		for (int i = numExamples * 5 / 6; i < numExamples; ++i) {
			TrainExamples.SetElement(0, i, randomReal(.3, .5));
			TrainExamples.SetElement(1, i, randomReal(1.0, -1.));
			TrainLabels.SetComponent(i, 6);
		}

		if (state == RenderState::Render3D) {
			for (int i = 0; i < numExamples * 1 / 6; ++i)
				TrainExamples.SetElement(2, i, randomReal(.4, 1.));
			for (int i = numExamples * 1 / 6; i < numExamples * 2 / 6; ++i)
				TrainExamples.SetElement(2, i, randomReal(.7, -1.));
			for (int i = numExamples * 2 / 6; i < numExamples * 3 / 6; ++i)
				TrainExamples.SetElement(2, i, randomReal(-0.6 / 1.3, 1.0 / 1.3));
			for (int i = numExamples * 3 / 6; i < numExamples * 4 / 6; ++i)
				TrainExamples.SetElement(2, i, randomReal(-0.5 / 1.3, -1.0 / 1.3));
			for (int i = numExamples * 4 / 6; i < numExamples * 5 / 6; ++i)
				TrainExamples.SetElement(2, i, randomReal(-.9, .9));
			for (int i = numExamples * 5 / 6; i < numExamples; ++i)
				TrainExamples.SetElement(2, i, randomReal(-1.0 / 1.5, -1.5 / 1.5));
		}
	}

	// Every vector ends with 1
	for (int i = 0; i < numExamples; ++i)
		TrainExamples.SetElement(numberOfDimensions, i, 1.0);
}

void SVM::MSVMAdapter::MakeImageDataPoints() {
	// const int numExamples = 10;
	// const int dimension = 180 * 300;
	// TrainExamples = Matrix<double>(dimension + 1, numExamples, true);
	// TestExamples = Matrix<double>(dimension + 1, numExamples, true);
	// TrainLabels = Vector<int>(numExamples);
	// double color = 0.;
	// for (int img = 1; img <= numExamples; img ++)
	// {
	// 	const png::image< png::rgb_pixel > image("../SVM/Data/" + to_string(img) + ".png");
	// 	for (int j = 0; j < 180; ++j)
	// 	{
	// 		for (int i = 0; i < 300; ++i)
	// 		{
	// 			color = image.get_pixel(j, i).red/255.;
	// 			TrainExamples.SetElement(j*300 + i, img-1, color);
	// 		}
	// 	}
	// }
	//
	// TrainLabels.SetComponent(0, 1);
	// TrainLabels.SetComponent(1, 1);
	// TrainLabels.SetComponent(2, 2);
	// TrainLabels.SetComponent(3, 2);
	// TrainLabels.SetComponent(4, 2);
	// TrainLabels.SetComponent(5, 2);
	// TrainLabels.SetComponent(6, 2);
	// TrainLabels.SetComponent(7, 1);
	// TrainLabels.SetComponent(8, 1);
	// TrainLabels.SetComponent(9, 1);
	//
	// // Every vector ends with 1
	// for (int i = 0; i < numExamples; ++i)
	// 	TrainExamples.SetElement(dimension, i, 1.0);
}

void SVM::MSVMAdapter::SetDataSet(int i) {
	currentDataSet = i;
}
