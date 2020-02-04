#include "Initialisation.h"

using namespace cv;
using namespace std;
using namespace Eigen;

VectorizationData initialisation(Image<Vec3b> I) {

	//INITIALIZE DATA
	Vec3b* C = new Vec3b(0, 0, 0);
	Image<Vec3b> image_copy;
	I.copyTo(image_copy);
	vector<Point2d>* vector_points = new vector<Point2d>;

	Image<Vec3b> colors = imread("../data/colors.jpg");

	VectorPoints data;
	data.image = I;
	data.image_copy = image_copy;
	data.C = C;
	data.colors = colors;
	data.vector_points = vector_points;

	resizeWindow(WINDOW_NAME, 600, 600);
	namedWindow(WINDOW_NAME, WINDOW_NORMAL);

	imshow(WINDOW_NAME, data.image);
	imshow(COLOR_NAME, data.colors);

	setMouseCallback(WINDOW_NAME, onMouse, &data);
	setMouseCallback(COLOR_NAME, onMouse_Color, &data);
	waitKey(0);
	destroyWindow(COLOR_NAME);

	cout << vector_points->size() << " points have been entered" << endl;
	Bezigon firstBezigon = Bezigon(*vector_points, *C, 0);
	firstBezigon.plot_curve(I, "Initialisation");

	VectorizationData initData = VectorizationData(firstBezigon, I);
	return initData;
};


void onMouse(int event, int x, int y, int flags, void* p) {
	if (event != EVENT_LBUTTONDOWN)
		return;
	Point m1(x, y);
	VectorPoints* data = (VectorPoints*)p;
	circle(data->image, m1, 2, Scalar(255, 0, 0), 1);
	resizeWindow(WINDOW_NAME, 600, 600);
	imshow(WINDOW_NAME, data->image);
	if (event == EVENT_LBUTTONDOWN)
	{
		Point end;
		if (data->vector_points->size() == 0) {
			cout << "Left button of the mouse is clicked - position (" << x << ", " << y << ")" << endl;
			data->vector_points->push_back(m1);
		}
		else {
			end = (*(data->vector_points)).back();
			if (end != m1) {
				cout << "Left button of the mouse is clicked - position (" << x << ", " << y << ")" << endl;
				data->vector_points->push_back(m1);
			}
			else { cout << "Please enter another point" << endl; }
		}
	}
}


void onMouse_Color(int event, int x, int y, int flags, void* p) {
	if (event != EVENT_LBUTTONDOWN)
		return;
	VectorPoints* data = (VectorPoints*)p;
	Point2d m1(x, y);

	circle(data->image_copy, m1, 4, Scalar(255, 0, 0), 2);
	imshow(COLOR_NAME, data->colors);
	*(data->C) = data->colors(x, y);
	cout << "Chosen color : " << data->colors(x, y) << endl;
}