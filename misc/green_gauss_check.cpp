#include <iostream>
#include <vector>
#include <cmath>

//---------------------------------------------------------
// Node class
//---------------------------------------------------------
class node {

  public:
    node() {};
    node(double xv, double yv) { x = std::vector<double>({xv,yv}); };
    double get_coord(const int i) const { return x[i]; };

  private:
    std::vector<double> x;
};

//---------------------------------------------------------
// Edge class
//---------------------------------------------------------
class edge {

  public:
    edge() {};
    edge(const node& n1, const node& n2);
    double get_midpoint(const int i) const { return xm[i]; };
    double get_length() const { return length; };

  private:
    std::vector<double> xm;
    double length;
};

// CTOR for edge
edge::edge(const node& n1, const node& n2) {

  // Computing the length of the edge
  length = sqrt(pow(n2.get_coord(0)-n1.get_coord(0),2) + pow(n2.get_coord(1)-n1.get_coord(1),2));

  // Computing the midpoint
  double xmid,ymid;
  xmid = 0.5*(n1.get_coord(0) + n2.get_coord(0));
  ymid = 0.5*(n1.get_coord(1) + n2.get_coord(1));
  xm = std::vector<double>({xmid,ymid});

}

//---------------------------------------------------------
// Element class
//---------------------------------------------------------
class element {

  public:
    element() {};
    element(std::vector<node>& n);
    std::vector<node> get_nodes() const { return nodes; };
    double get_centroid(const unsigned int i) const { return xc[i]; };
    void set_f(const double fv) { f = fv; };
    double get_f() const { return f; };

  private:
    std::vector<std::vector<double> > nhat;
    std::vector<node> nodes;
    double area;
    double f;  // function value defined at the cell center
    std::vector<double> xc;

};

// CTOR for element
element::element(std::vector<node>& n) {

  // Copying nodes
  nodes = n;
  xc.resize(2);

  // Computing area
  double x1,x2,x3,x4;
  double y1,y2,y3,y4;
  x1 = nodes[0].get_coord(0);
  x2 = nodes[1].get_coord(0);
  x3 = nodes[2].get_coord(0);
  x4 = nodes[3].get_coord(0);
  y1 = nodes[0].get_coord(1);
  y2 = nodes[1].get_coord(1);
  y3 = nodes[2].get_coord(1);
  y4 = nodes[3].get_coord(1);
  area = 0.5*((x1-x3)*(y2-y4)+(x4-x2)*(y1-y3));

  // Computing the centroid of the element
  xc[0] = 0.25*(x1+x2+x3+x4);
  xc[1] = 0.25*(y1+y2+y3+y4);

  // Computing face vectors
  nhat.resize(4);
  nhat[0] = std::vector<double>({y2-y1,x1-x2});
  nhat[1] = std::vector<double>({y3-y2,x2-x3});
  nhat[2] = std::vector<double>({y4-y3,x3-x4});
  nhat[3] = std::vector<double>({y1-y4,x4-x1});

  // Normalizing so that they are normal vectors
  double mag;
  for (int i=0; i<nhat.size(); ++i) {
    mag = sqrt(nhat[i][0]*nhat[i][0] + nhat[i][1]*nhat[i][1]);
    nhat[i][0] /= mag;
    nhat[i][1] /= mag;
  }

}

//---------------------------------------------------------
// Mesh class
//---------------------------------------------------------
class mesh {

  public:
    mesh() {};
    mesh(const unsigned int iMax, const unsigned int jMax, const double xmin, const double xmax, const double ymin, const double ymax);
    unsigned int get_imax() const { return imax; };
    unsigned int get_jmax() const { return jmax; };
    unsigned int get_nelemi() const { return nelemi; };
    unsigned int get_nelemj() const { return nelemj; };
    std::vector<std::vector<element> > get_elements() const { return elements; };
    void set_elements(const std::vector<std::vector<element> >& e) { elements = e; };

  private:
    unsigned int imax;
    unsigned int jmax;
    unsigned int nelemi;
    unsigned int nelemj;
    std::vector<std::vector<element> > elements;
    std::vector<std::vector<edge> > edges_h;    // horizontal edges
    std::vector<std::vector<edge> > edges_v;    // vertical edges

};

// CTOR for mesh
mesh::mesh(const unsigned int iMax, const unsigned int jMax, const double xmin, const double xmax, const double ymin, const double ymax) {

  // Declaring some variables
  imax = iMax;
  jmax = jMax;
  nelemi = imax-1;
  nelemj = jmax-1;
  std::vector<std::vector<double> > x(imax);
  std::vector<std::vector<double> > y(imax);
  std::vector<double> xtmp(jmax);
  std::vector<double> ytmp(jmax);

  // Figuring out grid spacing
  double dx,dy;
  dx = (xmax-xmin)/double(imax-1);
  dy = (ymax-ymin)/double(jmax-1);

  // Creating grid points
  for (unsigned int i=0; i<imax; ++i) {
    for (unsigned int j=0; j<jmax; ++j) {
      xtmp[j] = dx*i + xmin;
      ytmp[j] = dy*j + ymin;
    }
    x[i] = xtmp;
    y[i] = ytmp;
  }

  // Creating nodes
  std::vector<std::vector<node> > nodes(imax);
  std::vector<node> nodetmp(jmax);
  for (unsigned int i=0; i<imax; ++i) {
    for (unsigned int j=0; j<jmax; ++j) {
      nodetmp[j] = node(x[i][j],y[i][j]);
    }
    nodes[i] = nodetmp;
  }

  // Creating elements
  elements.resize(nelemi);
  std::vector<element> elemtmp(nelemj);
  std::vector<node> elemnodes(4);  // Assuming quads
  for (unsigned int i=0; i<nelemi; ++i) {
    for (unsigned int j=0; j<nelemj; ++j) {
      elemnodes[0] = nodes[i][j];
      elemnodes[1] = nodes[i+1][j];
      elemnodes[2] = nodes[i+1][j+1];
      elemnodes[3] = nodes[i][j+1];
      elemtmp[j] = element(elemnodes);
    }
    elements[i] = elemtmp;
  }

  // Creating horizontal edges
  edges_h.resize(imax);
  std::vector<edge> edgetmp(jmax);
  for (int i=0; i<nelemi; ++i) {
    for (int j=0; j<jmax; ++j) {
      edgetmp[j] = edge(nodes[i][j],nodes[i+1][j]);
    }
    edges_h[i] = edgetmp;
  }

  // Creating vertical edges
  edges_v.resize(imax);
  for (int i=0; i<imax; ++i) {
    for (int j=0; j<nelemj; ++j) {
      edgetmp[j] = edge(nodes[i][j],nodes[i][j+1]);
    }
    edges_v[i] = edgetmp;
  }

}

//---------------------------------------------------------
// Function for computing the gradient using the
// Green-Gauss approach
//---------------------------------------------------------
void compute_gradient_green_gauss(const mesh g, std::vector<std::vector<double> >& dfdx, std::vector<std::vector<double> >& dudy) {

  // Making sure 

}

int main() {

  // Creating a simple mesh
  mesh grid(50,50,0.0,1.0,0.0,1.0);
  std::cout << "Done generating mesh." << std::endl;

  // Defining an analytic function
  auto func = [] (double xv, double yv) { return sin(2.0*M_PI*xv)*sin(2.0*M_PI*yv); };

  // Filling in the values of the analytic function at the cell centers
  std::vector<std::vector<element> > elems = grid.get_elements();
  for (int i=0; i<grid.get_nelemi(); ++i) {
    for (int j=0; j<grid.get_nelemj(); ++j) {
      elems[i][j].set_f(func(elems[i][j].get_centroid(0),elems[i][j].get_centroid(1)));
    }
  }

  // Updating elements
  grid.set_elements(elems);

  // Computing the gradient using the Green-Gauss theorem

  return 0;

}
