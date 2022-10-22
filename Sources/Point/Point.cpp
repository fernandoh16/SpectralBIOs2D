#include <cmath> 
#include <cassert>
#include <iostream>
#include "Point.hpp"

Point::Point() {
	x[0] = 0; 
	x[1] = 0; 
	x[2] = 0;
	Tag  = 0;
}

Point::~Point() {

}

Point::Point(const Point & P) { 
	x[0] = P.x[0]; 
	x[1] = P.x[1]; 
	x[2] = P.x[2];
	Tag  = P.Tag;
}

Point::Point(double X, double Y, double Z, int otherTag)  { 
	x[0] = X; 
	x[1] = Y; 
	x[2] = Z;
	Tag = 0;
}

double & Point::operator[](int i) {
	assert(i >=0 ); 
	assert(i < 3);
	return x[i];
}

Point & Point::operator=(const Point& otherPoint) {
	x[0] = otherPoint.x[0]; 
	x[1] = otherPoint.x[1]; 
	x[2] = otherPoint.x[2]; 
	Tag  = otherPoint.Tag;
	return *this;	
} 

Point Point::operator+() const {
	Point v; 
	v[0] = x[0]; v[1] = x[1]; v[2] = x[2]; 
	return v;
}

Point Point::operator-() const {
	Point v; 
	v[0] = -x[0]; 
	v[1] = -x[1]; 
	v[2] = -x[2]; 
	return v;
}

Point Point::operator+(const Point & P) const {
	Point v; 
	v[0] = x[0] + P.x[0]; 
	v[1] = x[1] + P.x[1]; 
	v[2] = x[2] + P.x[2];
	return v;
}

Point Point::operator-(const Point & P) const {
	Point v; 
	v[0] = x[0] - P.x[0]; 
	v[1] = x[1] - P.x[1]; 
	v[2] = x[2] - P.x[2]; 
	return v;
}

Point Point::operator*(double a) const {
	Point v; 
	v[0] = a*x[0]; 
	v[1] = a*x[1]; 
	v[2] = a*x[2];
	return v;
}

double length(const Point & P) {
	return sqrt(pow(P.x[0],2)+pow(P.x[1],2)+pow(P.x[2],2));
}

double length(const Point & P1, const Point & P2) {
	return sqrt(pow(P1.x[0]-P2.x[0],2)+pow(P1.x[1]-P2.x[1],2)+pow(P1.x[2]-P2.x[2],2));
}

Point Point::operator*(const Point & P1) const{
	Point v;
	v[0] = -P1.x[1]*x[2] + P1.x[2]*x[1]; 
	v[1] = -P1.x[2]*x[0] + P1.x[0]*x[2]; 
	v[2] = -P1.x[0]*x[1] + P1.x[1]*x[0]; 
	return v;
}

double Point::operator&(const Point & P1) const{	
	return P1.x[0]*x[0]+P1.x[1]*x[1]+P1.x[2]*x[2]; 
}

int Point::GetTag() {
	return Tag;
}

void Point::SetTag(int Tag_) {
	Tag = Tag_;
}