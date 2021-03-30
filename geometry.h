#ifndef GEOMETRY_GEOMETRY_H
#define GEOMETRY_GEOMETRY_H
#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>

const double EPS = 0.000000001;

bool is_equal (double a, double b) {
    return std::fabs(a-b) < EPS;
}
class Line;

struct Point {
    double x;
    double y;
    Point (double a, double b) {
        x = a;
        y = b;
    }
    Point() {}
    Point (const Point& p) {
        x = p.x;
        y = p.y;
    }
    Point& operator=(const Point& p) {
        x = p.x;
        y = p.y;
        return *this;
    }
    double getDistance(Point p) const {
        return sqrt((x-p.x)*(x-p.x) + (y-p.y)*(y-p.y));
    }
    Point getMiddle(Point p) const {
        return Point((p.x+x)/2, (p.y+y)/2);
    }

    Point getScaled(Point O, double k) const;

    void rotate (Point p, double angle) {
        Point v (x - p.x, y - p.y);
        double x1 = v.x;
        double y1 = v.y;
        v.x = v.x*cos(angle) - v.y*sin(angle);
        v.y = y1*cos(angle) + x1*sin(angle);
        x = p.x + v.x;
        y = p.y + v.y;
    }
    Point gerReflexed (Line) const;
};

bool operator== (const Point& a, const Point& b) {
    return (is_equal(a.x, b.x) && is_equal(a.y, b.y));
}

bool operator!= (const Point& a, const Point& b) {
    return !(a==b);
}

class Line { //ax+by+c = 0
public:
    double a;
    double b;
    double c;
    Line (Point a, Point b) {
        if (is_equal(a.x, b.x)) {
            this->a  = 1;
            this->b = 0;
            c = -a.x;
            return;
        }
        this->a = (b.y - a.y)/(b.x-a.x);
        this->c = b.y-this->a*b.x;
        this->b = -1;
    }
    Line (double coef, double mv) {
        a = coef;
        c = mv;
        b = -1;
    }
    Line (Point p, double coef) {
        a = coef;
        c = p.y - a*p.x;
        b = -1;
    }
    Line (double aa, double bb, double cc) {
        a = aa;
        b = bb;
        c = cc;
    }
    Line() = default;

    Point getIntersection(Line l) const;

    bool intersects(Line l) const;

    Line (const Line& l) {
        a = l.a;
        b = l.b;
        c = l.c;
    }

    Line& operator= (const Line& l) {
        a = l.a;
        b = l.b;
        c = l.c;
        return *this;
    }
    bool contains(Point) const;
};

Point Point::gerReflexed(Line axis) const {
    double t = (-axis.a*x - axis.b*y - axis.c)/ (axis.a * axis.a + axis.b*axis.b);
    return Point(2*axis.a*t + x, 2*axis.b*t+y);
}

bool operator== (const Line& a, const Line& b) {
    if (!is_equal(a.a, 0) && !is_equal(b.a, 0)) return (is_equal(a.b/a.a,b.b/b.a) && is_equal(a.c/a.a, b.c/b.a));
    else if (is_equal(a.a, 0) && is_equal(b.a, 0)) return is_equal(a.c/a.b, b.c/b.b);
    else return false;
}
bool operator!= (const Line& a, const Line& b) {
    return !(a==b);
}


class Shape {
public:
    virtual double perimeter() const = 0;
    virtual double area() const = 0;
    virtual bool operator==(const Shape& another)const  = 0;
    virtual bool isCongruentTo(const Shape& another) const = 0;
    virtual bool isSimilarTo(const Shape& another) const = 0;
    virtual bool containsPoint(Point point) const = 0;
    virtual void rotate(Point center, double angle) =  0;
    virtual void reflex(Point center) = 0;
    virtual void reflex(Line axis) = 0;
    virtual void scale(Point center, double coefficient) = 0;
    virtual ~Shape() = default;
};

Point Point::getScaled (Point O, double k) const {
    double Ox = x - O.x;
    double Oy = y - O.y;
    return Point(O.x + Ox*k, O.y + Oy*k);
}

class Polygon : public Shape {
protected:
    std::vector<Point> vertices;
public:
    size_t verticesCount() const {
        return vertices.size();
    }
    const std::vector<Point>& getVertices () const {
        return vertices;
    }
    bool isConvex() const;
    explicit Polygon (std::vector<Point>& p) {
        vertices = p;
    }
    Polygon() {}
    void addVertices() {}
    template <class Head, class... Tail>
    void addVertices (const Head& head, const Tail&...tail) {
        vertices.push_back(head);
        addVertices(tail...);
    }

    template <class Head, class... Tail>
    Polygon (const Head& head, const Tail&...tail) {
        addVertices(head, tail...);
    }

    void scale(Point center, double coefficient);
    void reflex(Point center) override {
        scale(center, -1);
    }

    double area() const;

    bool containsPoint(Point p) const;

    double perimeter() const override;

    bool operator==(const Shape& other) const;

    bool isSimilarTo(const Shape &another) const;

    std::vector<double> getAngles() const;

    bool isCongruentTo(const Shape &another) const;

    std::vector<double> getSides() const;

    void rotate(Point center, double angle);

    void reflex(Line axis);
    bool operator!= (const Shape& other) const;

    double getCoef(const Shape &another) const;
};

bool Line::contains (Point p) const {
    return is_equal(a*p.x+b*p.y+c, 0);
}

double Polygon::area() const {
    double s = 0;
    for (size_t i = 0; i< vertices.size() - 1; ++i) {
        s += vertices[i].x * vertices[i+1].y;
        s -= vertices[i].y * vertices[i+1].x;
    }
    s+= vertices[0].y * vertices[verticesCount() -1].x;
    s-=vertices[0].x * vertices[verticesCount() -1].y;
    s = std::fabs(s)/2;
    return s;
}

bool between(double a, double b, double c) {
    return (a<=b && b<= c) || (c<=b && b<= a);
}

bool Polygon::containsPoint(Point p) const {
    Line l(p, Point(p.x+1, p.y));
    int intersections = 0;
    for (size_t i = 0; i< vertices.size(); ++i) {
        Line edge(vertices[i], vertices[(i+1)%vertices.size()]);
        if (edge.intersects(l)) {
            Point inter = edge.getIntersection(l);
            if (inter == vertices[i] || inter == vertices[(i+1)%verticesCount()]) {
                bool intr = ((vertices[i].y < p.y && vertices[(i + 1) % vertices.size()].y >= p.y)
                            || (vertices[i].y >= p.y && vertices[(i + 1) % vertices.size()].y < p.y)) && inter.x >= p.x;
                if (intr) ++intersections;
            } else if (between(vertices[i].y, inter.y, vertices[(i+1)%verticesCount()].y) && inter.x>=p.x)
                    ++intersections;
        }
        if (edge.contains(p)) {
            if (between(vertices[i].x, p.x, vertices[(i+1)%verticesCount()].x)) return true;
        }
    }
    return (intersections%2!=0);
}

bool Polygon::isConvex() const {
    assert(verticesCount() >= 3);
    Point ab(vertices[1].x - vertices[0].x ,vertices[1].y - vertices[0].y);
    Point bc(vertices[2].x - vertices[1].x ,vertices[2].y - vertices[1].y);
    double d = ab.x*bc.y - ab.y*bc.x;
    assert(d!=0); // if d==0 the polygon has three point lying on the same line
    bool sign = d>0;
    for (size_t i = 1; i<vertices.size(); ++i) {
        ab = Point(vertices[(i+1)%verticesCount()].x - vertices[i].x ,
                vertices[(i+1)%verticesCount()].y - vertices[i].y);
        bc = Point(vertices[(i+2)%verticesCount()].x - vertices[(i+1)%verticesCount()].x ,
                vertices[(i+2)%verticesCount()].y - vertices[(i+1)%verticesCount()].y);
        d = ab.x*bc.y - ab.y*bc.x;
        assert(d!=0);
        if ((d>0) != sign) return false;
    }
    return true;
}

std::vector<double> Polygon::getSides () const {
    std::vector<double> lens;
    for (size_t i = 0; i< vertices.size(); ++i)
        lens.push_back(vertices[i].getDistance(vertices[(i+1)%vertices.size()]));
    return lens;
}
double Polygon::perimeter() const{
    double p = 0;
    auto sides = getSides();
    for (auto side: sides) p+=side;
    return p;
}


bool Polygon::operator==(const Shape& other) const {
    const Polygon* oth = dynamic_cast<const Polygon*>(&other);
    if (!oth) return false;
    if (oth->verticesCount() != verticesCount()) return false;
    int ind = -1;
    for (size_t i = 0; i<verticesCount(); ++i)
        if (vertices[i] == oth->vertices[0]) {
            ind  = i;
            break;
        }
    if (ind == -1) return false;
    int step = 0;
    if (vertices[(ind+1)%verticesCount()] == oth->vertices[1]) step = 1;
    else if (vertices[(verticesCount() + ind - 1)%verticesCount()] == oth->vertices[1]) step = -1;
    if (step == 0) return false;
    int new_ind = ind + step;
    for (size_t i = 2; i<verticesCount(); ++i) {
        new_ind += step;
        if (new_ind<0) new_ind += verticesCount();
        if (static_cast<size_t>(new_ind) >= verticesCount()) new_ind -= verticesCount();
        if (oth->vertices[i] != vertices[new_ind]) return false;
    }
    return true;
};

bool Polygon::operator!= (const Shape& other) const {
    return !(operator==(other));
}

std::vector<double> Polygon::getAngles () const {
    int n = verticesCount();
    std::vector<double> angles;
    for (size_t i = 0; i<verticesCount(); ++i) {
        Point b = vertices[i];
        Point c = vertices[(i+1)%n];
        Point a = vertices[(i-1+n)%n];
        Point v1 (a.x - b.x, a.y - b.y), v2(c.x - b.x, c.y - b.y);
        angles.push_back(((v1.x*v2.x + v1.y*v2.y)/(a.getDistance(b)*b.getDistance(c))));
    }
    return angles;
}

double Polygon::getCoef (const Shape& another) const {
    auto other = dynamic_cast<const Polygon*>(&another);
    if (!other) return 0;
    if (verticesCount() != other->verticesCount()) return 0;
    auto angles1 = getAngles(), angles2 = other->getAngles();
    auto sides1 = getSides(), sides2 = other->getSides();
    std::vector<int> indexes;
    for (int i = 0; i< verticesCount(); ++i)
        if (is_equal(angles1[i], angles2[0])) indexes.push_back(i);
    double k = 0;
    for (int ind: indexes) {
        int step = 0;
        if (is_equal(angles1[(ind+1)%verticesCount()], angles2[1])) step = 1;
        else if (is_equal(angles1[(verticesCount() + ind - 1)%verticesCount()], angles2[1])) step = -1;
        if (step == 0) continue;
        int new_ind = ind;
        if (step == -1) k = sides1[(new_ind-1+verticesCount())%verticesCount()]/ sides2[0];
        else k = (sides1[new_ind] / sides2[0]);
        int edge_start_ind;
        for (size_t i = 1; i<verticesCount(); ++i) {
            if (step == 1) edge_start_ind = new_ind;
            new_ind += step;
            if (new_ind>= static_cast<int>(verticesCount())) new_ind -= verticesCount();
            if (new_ind<0) new_ind += verticesCount();
            if (step == -1) edge_start_ind = new_ind;
            if (!is_equal(angles2[i], angles1[new_ind])) k = 0;
            if (!is_equal(k, sides1[edge_start_ind] / sides2[i-1])) k = 0;
            if (is_equal(k, 0)) break;
        }
        if (!is_equal(k, 0)) return k;
    }
    return k;
}

bool Polygon::isSimilarTo(const Shape &another) const {
    return (!is_equal(getCoef(another), 0));
}

bool Polygon::isCongruentTo(const Shape &another) const {
    return is_equal(getCoef(another), 1);
}

void Polygon::scale(Point center, double coefficient) {
    for (size_t i = 0; i<vertices.size(); ++i)
        vertices[i] = vertices[i].getScaled(center, coefficient);
}

void Polygon::rotate(Point center, double angle) {
    angle *= (M_PI/180);
    for (size_t i = 0; i<verticesCount(); ++i)
        vertices[i].rotate(center, angle);
}

void Polygon::reflex(Line axis) {
    for (size_t i = 0; i<vertices.size(); ++i)
        vertices[i] = vertices[i].gerReflexed(axis);
}

class Ellipse : public Shape {
protected:
    double b;
    Point f1;
    Point f2;
    double a;
public:
    std::pair<Point,Point> focuses() const {
        return std::make_pair(f1, f2);
    }
    Line getDirectice (Point f1) const {
        double e = sqrt(a*a - b*b) / a;
        Point v (f1.x - center().x, f1.y - center().y);
        v.x *= (1/e);
        v.y *= (1/e);
        Point v2 = v;
        v.x = center().x + v.x;
        v.y = center().y + v.y;
        Point p (1,-v2.x/v2.y);
        Point m1 (center().x + p.x + v.x, center().y + p.y + v.y);
        return Line(v, m1);
    }
    std::pair<Line, Line> directrices() const {
        return std::make_pair(getDirectice(f1), getDirectice(f2));
    }
    double eccentricity() const {
        return sqrt(a*a - b*b) / a;
    }

    virtual Point center () const {
        return f1.getMiddle(f2);
    }
    Ellipse (Point F1, Point F2, double dist) {
        f1 = F1;
        f2 = F2;
        a = dist/2;
        double c = center().getDistance(f1);
        b  = sqrt(a*a - c*c);
    }
    Ellipse() = default;
    void scale(Point center, double coefficient) override;
    void reflex(Point center) override {
        scale(center, -1);
    }
    bool operator== (const Shape& other) const override ;

    bool isSimilarTo(const Shape &another) const override;

    bool isCongruentTo(const Shape &another) const override;

    bool containsPoint(Point point) const override;

    double area() const override;

    double perimeter() const override;

    void rotate(Point center, double angle);

    void reflex(Line axis);
    bool operator!= (const Shape&) const;
};

bool Ellipse::operator== (const Shape& other) const {
    const Ellipse* oth = dynamic_cast<const Ellipse*>(&other);
    if (!oth) return false;
    return (is_equal(a, oth->a) && f1 == oth->f1 && f2 == oth->f2);
}

bool Ellipse::operator!= (const Shape& other) const {
    return !(operator==(other));
}

bool Ellipse::isCongruentTo(const Shape &another) const {
    const Ellipse* oth = dynamic_cast<const Ellipse*>(&another);
    if (!oth) return false;
    return (is_equal(a, oth->a) && is_equal(b, oth->b));
}

bool Ellipse::isSimilarTo(const Shape &another) const {
    const Ellipse* oth = dynamic_cast<const Ellipse*>(&another);
    if (!oth) return false;
    return is_equal((b/a), (oth->b/oth->a));
}

bool Ellipse::containsPoint(Point point) const { //not sure if it is correct
    return (point.getDistance(f1)+point.getDistance(f2) <= 2*a + EPS);
}

void Ellipse::scale(Point center, double coefficient) {
    f1 = f1.getScaled(center, coefficient);
    f2 = f2.getScaled(center, coefficient);
    a = std::fabs(coefficient*a);
    b = std::fabs(coefficient*b);
}

double Ellipse::area() const {
    return M_PI*a*b;
}

double Ellipse::perimeter() const {
    return M_PI*(3*(a+b)- sqrt((3*a+b)*(a+3*b)));
}

void Ellipse::rotate(Point center, double angle) {
    angle *= (M_PI/180);
    f1.rotate(center, angle);
    f2.rotate(center, angle);
}

void Ellipse::reflex(Line axis) {
    f1 = f1.gerReflexed(axis);
    f2 = f2.gerReflexed(axis);
}

Point Line::getIntersection(Line l) const {
    Point p;
    assert(a*l.b != l.a*b);
    p.x = (l.c*b - c*l.b)/(a*l.b - l.a* b);
    if (l.b !=0) p.y = (-l.c - l.a*p.x)/l.b;
    else p.y = (c -a*p.x)/b;
    return p;
}

bool Line::intersects(Line l) const {
    return a*l.b != l.a*b;
}


class Circle : public Ellipse {
    Point centre;
public:
    double radius () const {
        return a;
    }
    Circle (Point c, double r) {
        centre = c;
        a = r;
        b = r;
    }
    double area () const override {
        return M_PI * a* a;
    }
    double perimeter () const override { return 2 * M_PI * a; }
    bool containsPoint (Point p) const override {
        return (p.getDistance(centre) <= radius());
    }

    void scale(Point center, double coefficient) override;
    void reflex(Point center) override {
        scale(center, -1);
    }
    Point center() const override { return centre; }
    bool operator== (const Shape& other) const override ;
    bool operator!= (const Shape& other) const;
    void rotate(Point center, double angle) override;

    void reflex(Line axis) override;

    Circle& operator= (const Circle& c) {
        centre = c.centre;
        a = c.a;
        b = c.b;
        return *this;
    }

    Circle (const Circle& c) {
        centre = c.centre;
        a = c.a;
        b = c.b;
    }
};

bool Circle::operator== (const Shape& another) const {
    const Circle* oth = dynamic_cast<const Circle*>(&another);
    if (!oth) return false;
    return (center() == oth->center() && is_equal(radius(), oth->radius()));
}
bool Circle::operator!= (const Shape& other) const {
    return !(operator==(other));
}

void Circle::scale(Point center, double coefficient) {
    centre = centre.getScaled(center, coefficient);
    a = std::fabs(coefficient*a);
    b = std::fabs(coefficient*b);
}

void Circle::rotate(Point center, double angle) {
    angle *= (M_PI/180);
    centre.rotate(center, angle);
}

void Circle::reflex(Line axis) {
    centre = centre.gerReflexed(axis);
}

class Rectangle : public Polygon {
public:
    std::pair<Line, Line> diagonals() const {
        return std::make_pair(Line(vertices[0], vertices[2]), Line(vertices[1], vertices[3]));
    }
    Point center() const {
        auto d = diagonals();
        return d.first.getIntersection(d.second);
    }
    Rectangle (Point a, Point b, Point c, Point d) : Polygon(a,b,c,d) {}
    Rectangle(Point a, Point b, double coef) {
        vertices.resize(4);
        if (coef>1) coef = 1/coef;
        double cos = (1-coef*coef)/(1+coef*coef);
        Point o = a.getMiddle(b);
        vertices[0] = a;
        a.rotate(o, -acos(cos));
        vertices[1] = a;
        vertices[2] = b;
        b.rotate(o, -acos(cos));
        vertices[3] = b;

    }
    explicit Rectangle(std::vector<Point>& points) : Polygon(points) {}

    double area() const override;

    double perimeter() const override;
};


double Rectangle::area() const {
    double a = vertices[0].getDistance(vertices[1]);
    double b = vertices[1].getDistance(vertices[2]);
    return a*b;
}

double Rectangle::perimeter() const {
    double a = vertices[0].getDistance(vertices[1]);
    double b = vertices[1].getDistance(vertices[2]);
    return 2*(a+b);
}

class Square : public Rectangle {
public:
    Square (Point a, Point b) : Rectangle(a,b,1.0) {}

    explicit Square(std::vector<Point>& points) : Rectangle(points) {}

    Square (Point a, Point b, Point c, Point d) : Rectangle(a,b,c,d) {}

    Circle inscribedCircle() const;

    double area() const override;

    double perimeter() const override;

    Circle circumscribedCircle() const;

};

Circle Square::inscribedCircle() const{
    return Circle(center(), vertices[0].getDistance(vertices[1])/2);
}
Circle Square::circumscribedCircle() const {
    return Circle(center(), vertices[0].getDistance(vertices[1]) * sqrt(2) / 2);
}

double Square::area() const {
    double a = vertices[0].getDistance(vertices[1]);
    return a*a;
}

double  Square::perimeter() const {
    return 4* vertices[0].getDistance(vertices[1]);
}

class Triangle : public Polygon {
public:
    explicit Triangle(std::vector<Point>& points) : Polygon(points) {}

    Triangle (Point a, Point b, Point c) : Polygon(a,b,c) {}

    Circle circumscribedCircle() const;

    Circle ninePointsCircle() const;

    Line EulerLine() const;

    Point orthocenter() const;

    Line getHeight(int base_index) const;

    Line getMidPerpendicular(Point a, Point b) const;

    Point centroid() const;

    Circle inscribedCircle() const;

    void print() const;
};

Line getBisector (Point a, Point b, Point c) { // a is the angle's vertex (angle is BAC)
    double coef = a.getDistance(b) / a.getDistance(c);
    Point k ((b.x+coef*c.x)/(1+coef), (b.y+coef*c.y)/(1+coef));
    return Line(a, k);
}

Line Triangle::getMidPerpendicular (Point a, Point b) const {
    Point m = a.getMiddle(b);
    Point v (b.x-a.x, b.y - a.y);
    Point p;
    if (v.y != 0) p = Point(1,-v.x/v.y);
    else p = Point(0,1);
    Point m1 (p.x + m.x, p.y + m.y);
    return Line(m, m1);
}
Circle Triangle::circumscribedCircle() const {
    auto sides = getSides();
    return Circle (getMidPerpendicular(vertices[0], vertices[1]).getIntersection(
            getMidPerpendicular(vertices[1], vertices[2])), sides[0]*sides[1]*sides[2] / (4*area()));
}

Circle Triangle::inscribedCircle () const {
    Point center = getBisector(vertices[0], vertices[1], vertices[2]).getIntersection(
            getBisector(vertices[1], vertices[0], vertices[2]));
    return Circle(center, area()*2/perimeter());
}


Line Triangle::getHeight (int base_index) const {
     return Line(vertices[base_index],
             vertices[base_index].gerReflexed(
                     Line(vertices[(base_index+1)%3], vertices[(base_index+2)%3])));
}

Circle Triangle::ninePointsCircle() const {
    std::vector<Point> middles(3);
    for (int i = 0; i<3; ++i) {
        middles[i] = vertices[i].getMiddle(vertices[(i+1)%3]);
    }
    return Triangle(middles).circumscribedCircle();
}

Line Triangle::EulerLine() const {
    return Line(circumscribedCircle().center(), orthocenter());
}

Point Triangle::orthocenter() const {
    return getHeight(0).getIntersection(getHeight(1));
}

Point Triangle::centroid () const {
    Line m1 (vertices[0], vertices[1].getMiddle(vertices[2]));
    Line m2 (vertices[1], vertices[0].getMiddle(vertices[2]));
    return m1.getIntersection(m2);
}


#endif //GEOMETRY_GEOMETRY_H