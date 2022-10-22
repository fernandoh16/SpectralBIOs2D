#ifndef POINTHEADERDEF 
#define POINTHEADERDEF 

class Point {

private:
	double x[3];
	int Tag;
	
public:
	Point();
	~Point();
	Point(const Point & P);
	Point(double X, double Y, double Z = 0, int otherTag = 0);
	double & operator[](int i); 
	Point & operator = (const Point & otherPoint); 
	Point operator+() const; 
	Point operator-() const; 
	Point operator+(const Point & P) const; 
	Point operator-(const Point & P) const; 
	Point operator*(double a) const; 
	friend double length(const Point & P); 
	friend double length(const Point & P1, const Point & P2); 
	Point operator*(const Point & P1) const; 
	double operator&(const Point & P1) const;
	int GetTag();
	void SetTag(int Tag_);
};
	double length(const Point & P); 
	double length(const Point & P1, const Point & P2); 
	 
#endif