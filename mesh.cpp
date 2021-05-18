#include "mesh.h"
#include "matrix.h"
#include <cstring>
#include <iostream>
#include <strstream>
#include <fstream>
#include <cmath>
#include <float.h>
using namespace std;
constexpr auto PI = 3.1415926535;   // pi;

/////////////////////////////////////////
// helping inline functions
// cot<p2p1, p2p3>
inline double Cot(const Vector3d & p1, const Vector3d & p2, const Vector3d & p3) {
	Vector3d v1 = p1 - p2;
	Vector3d v2 = p3 - p2;

	v1 /= v1.L2Norm();
	v2 /= v2.L2Norm();
	double tmp = v1.Dot(v2);
	return 1.0 / tan(acos(tmp));
}

// Area(p1, p2, p3): cross<p1p2, p1p3>
inline double Area(const Vector3d & p1, const Vector3d & p2, const Vector3d & p3) {
	Vector3d v1 = p2 - p1;
	Vector3d v2 = p3 - p1;
	return v1.Cross(v2).L2Norm() / 2.0;
}


/////////////////////////////////////////
// implementation of OneRingHEdge class
OneRingHEdge::OneRingHEdge(const Vertex * v) {
	if (v == NULL) start = next = NULL;
	else start = next = v->HalfEdge();
}

HEdge * OneRingHEdge::NextHEdge() {
	HEdge *ret = next;
	if (next && next->Prev()->Twin() != start)
		next = next->Prev()->Twin();
	else
		next = NULL;
	return ret;
}

/////////////////////////////////////////
// implementation of Mesh class
//
// function AddFace
// it's only for loading obj model, you do not need to understand it
void Mesh::AddFace(int v1, int v2, int v3) {
	int i;
	HEdge *he[3], *bhe[3];
	Vertex *v[3];
	Face *f;

	// obtain objects
	for (i=0; i<3; i++) he[i] = new HEdge();
	for (i=0; i<3; i++) bhe[i] = new HEdge(true);
	v[0] = vList[v1];
	v[1] = vList[v2];
	v[2] = vList[v3];
	f = new Face();

	// connect prev-next pointers
	SetPrevNext(he[0], he[1]);
	SetPrevNext(he[1], he[2]);
	SetPrevNext(he[2], he[0]);
	SetPrevNext(bhe[0], bhe[1]);
	SetPrevNext(bhe[1], bhe[2]);
	SetPrevNext(bhe[2], bhe[0]);

	// connect twin pointers
	SetTwin(he[0], bhe[0]);
	SetTwin(he[1], bhe[2]);
	SetTwin(he[2], bhe[1]);

	// connect start pointers for bhe
	bhe[0]->SetStart(v[1]);
	bhe[1]->SetStart(v[0]);
	bhe[2]->SetStart(v[2]);
	for (i=0; i<3; i++) he[i]->SetStart(v[i]);

	// connect start pointers
	// connect face-hedge pointers
	for (i=0; i<3; i++) {
		v[i]->SetHalfEdge(he[i]);
		v[i]->adjHEdges.push_back(he[i]);
		SetFace(f, he[i]);
	}
	v[0]->adjHEdges.push_back(bhe[1]);
	v[1]->adjHEdges.push_back(bhe[0]);
	v[2]->adjHEdges.push_back(bhe[2]);

	// mearge boundary if in need
	for (i=0; i<3; i++) {
		Vertex *start = bhe[i]->Start();
		Vertex *end   = bhe[i]->End();
		for (size_t j=0; j<end->adjHEdges.size(); j++) {
			HEdge *curr = end->adjHEdges[j];
			if (curr->IsBoundary() && curr->End()==start) {
				SetPrevNext(bhe[i]->Prev(), curr->Next());
				SetPrevNext(curr->Prev(), bhe[i]->Next());
				SetTwin(bhe[i]->Twin(), curr->Twin());
				bhe[i]->SetStart(NULL);	// mark as unused
				curr->SetStart(NULL);	// mark as unused
				break;
			}
		}
	}

	// finally add hedges and faces to list
	for (i=0; i<3; i++) heList.push_back(he[i]);
	for (i=0; i<3; i++) bheList.push_back(bhe[i]);
	fList.push_back(f);
}

// function LoadObjFile
// it's only for loading obj model, you do not need to understand it
bool Mesh::LoadObjFile(const char *filename) {
	if (filename==NULL || strlen(filename)==0) return false;
	ifstream ifs(filename);
	if (ifs.fail()) return false;

	Clear();

	char buf[1024], type[1024];
	do {
		ifs.getline(buf, 1024);
		istrstream iss(buf);
		iss >> type;

		// vertex
		if (strcmp(type, "v") == 0) {
			double x, y, z;
			iss >> x >> y >> z;
            AddVertex(new Vertex(x,y,z));
		}
		// face
		else if (strcmp(type, "f") == 0) {
			int index[3];
			iss >> index[0] >> index[1] >> index[2];
			AddFace(index[0]-1, index[1]-1, index[2]-1);
		}
	} while (!ifs.eof());
	ifs.close();

	size_t i;
	Vector3d box = this->MaxCoord() - this->MinCoord();
	for (i=0; i<vList.size(); i++) vList[i]->SetPosition(vList[i]->Position() / box.X());

	Vector3d tot;
	for (i=0; i<vList.size(); i++) tot += vList[i]->Position();
	Vector3d avg = tot / vList.size();
	for (i=0; i<vList.size(); i++) vList[i]->SetPosition(vList[i]->Position() - avg);

	HEdgeList list;
	for (i=0; i<bheList.size(); i++)
		if (bheList[i]->Start()) list.push_back(bheList[i]);
	bheList = list;

	for (i=0; i<vList.size(); i++) 
	{
		vList[i]->adjHEdges.clear();
		vList[i]->SetIndex(static_cast<int> (i));
		vList[i]->SetFlag(0);
	}

	return true;
}

/////////// helper functions:///// 
// traverse connected edges
void visitBounds(HEdge *e) {
	while (!e->Next()->Flag() && e->Next()->IsBoundary()) {
		e->SetFlag(true);
		e = e->Next();
	}
}
//walk along all boundary edges
int Mesh::CountBoundaryLoops()
{
	int bounds = bheList.size(), loops = 0;
	if (bounds == 0) return 0;
	for (int i = 0; i < bounds; i++) {
		HEdge *edge = bheList[i];
		if (!edge->Flag()) {
			loops++;
			visitBounds(edge);
		}
	}
	// flip all flags back to flase
	for (int i = 0; i < bounds; i++) {
		HEdge* edge = bheList[i];
		edge->SetFlag(0);
	}
	return loops;
}

// dfs
void dfs(HEdge* e) {
	HEdgeList neighbors;
	neighbors.push_back(e);

	while (!neighbors.empty())
	{
		e = neighbors.back();
		neighbors.pop_back();
		e->SetFlag(1);

		if (e->Next() != NULL) {
			if (!(e->Next()->Flag())) {
				neighbors.push_back(e->Next());
			}
		}
		if (e->Twin() != NULL) {
			if (!(e->Twin()->Flag())) {
				neighbors.push_back(e->Twin());
			}
		}
	}
}
// dfs, bfs traverse to find how many independent polygons
int Mesh::CountConnectedComponents()
{
	size_t edges = heList.size(), cc = 0;

	if (edges <= 12) return 1;
	for (size_t i = 0; i < edges; i++) {
		HEdge *e = heList[i];
		if (!e->Flag()) {
			// cout << "visiting hedge: " << e << endl;
			cc++;
			dfs(e);
		}
	}
	// flip all flags back to flase
	for (size_t i = 0; i < edges; i++) {
		HEdge* e = heList[i];
		e->SetFlag(false);
	}
	return cc;
}
////////// helper function ends/////

void Mesh::DisplayMeshInfo()
{
	/*************************/
	/* insert your code here */
	/*************************/	
	int e = heList.size() / 2 + bheList.size(), v = vList.size(), f = fList.size();
	cout << "length of Hedge: " << heList.size() <<
		"\nlength of bHedge: "	<< bheList.size() <<
		"\nEdge number: "		<< e <<
		"\nVertex Number: "		<< v <<
		"\nFace number: "		<< f << endl;

	// Euler characteristic
	// calculate b,c
	int b = CountBoundaryLoops();
	int c = CountConnectedComponents();
	cout << "number of boundary loop: " << b <<
		"\nnumber of connected components: " << c << endl;

	// num of genus
	// v - e + f = 2(c - g) - b 	
	int g = c - (v - e + f + b) / 2;
	cout << "number of genus: "		<< g << endl;
}


// -------------------------------------------------------
// DO NOT TOUCH THE FOLLOWING FOR NOW
// -------------------------------------------------------
void Mesh::ComputeVertexNormals() 
{
	/*************************/
	/* insert your code here */
	/*************************/
	VertexList l = vList;
	for (size_t i = 0; i < l.size(); i++) {
		Vertex* p = l[i];
		int k = p->Valence(), j=0;
		OneRingVertex ring(p);
		Vertex* p_j = NULL;

		if (!p->IsBoundary()) {
			// interior vector
			Vector3d t1 = Vector3d(0, 0, 0), t2 = Vector3d(0, 0, 0);
			double theta = 2 * PI / k;

			while(p_j=ring.NextVertex()){	
				t1 += (p_j->Position()) * cos(j * theta);
				t2 += (p_j->Position()) * sin(j * theta);
				j++;
			}
			
			// update normal with the crross product
			Vector3d normal = t1.Cross(t2); 
			normal /= normal.L2Norm();
			p->SetNormal(normal);
		}else {
			// boundary vector 
			Vector3d t_along, t_across;
			double theta = PI / ((double)k - 1);
			while (p_j = ring.NextVertex()){
				if (j == 0) {
					// p_0
					t_along = p_j->Position();
					if (k == 2)	t_across = p_j->Position();
					if (k >= 4)	t_across = (p_j->Position()) * sin(theta);
				}
				if (j == 1) {
					// p_1
					if (k == 2) t_across += p_j->Position();
					if (k == 3) t_across = p_j->Position();
				}
				if (j == k - 1) {
					// p_{k-1}
					t_along -= p_j->Position();
					if (k >= 4) t_across += (p_j->Position()) * sin(theta);
				}
				if (j >= 1 && j <= k - 2 && k >= 4)
					t_across += (p_j->Position()) * sin(j * theta) * (2. * cos(theta) - 2.);

				j++;
			}
			if (k == 2) t_across -= (p->Position()) * 2;
			if (k == 3) t_across -= p->Position();

			// set cross product as the new normal
			Vector3d normal = t_along.Cross(t_across);
			normal /= normal.L2Norm();
			p->SetNormal(normal);
		}
	}
}

void Mesh::UmbrellaSmooth() 
{
	/*************************/
	/* insert your code here */
	/*************************/
	size_t v = vList.size();
	Matrix *factor = new Matrix(v, v);	// factor = I+lamba*L, L: cotangent-weight Matrix
	double* vx = new double[v], * vy = new double[v], * vz = new double[v];// old coord matrix
	double* vx1 = new double[v](), * vy1 = new double[v](), * vz1 = new double[v](), lamda = 0.9;// new coord matrix

	// update matrix
	for (size_t i = 0; i < v; i++) {// i-th row
		Vertex* p_i = vList[i], *p_j;
		int k = p_i->Valence();	// j-th column
		if (p_i->IsBoundary() || k < 3) continue;
		HEdge* start = p_i->HalfEdge();
		HEdge* next = start;
		Vector3d p_j_l, p_j_r;
		double wj_sum=0;
		if(i%500==0)cout << "index=" << p_i->Index() << " i=" << i << endl;
		// init old coord matrix
		vx[i] = p_i->Position().X();
		vy[i] = p_i->Position().Y();
		vz[i] = p_i->Position().Z();
		factor->AddElement(i, i, 1 - lamda);

		// calculate wij for each vertex
		do{
			p_j = next->Next()->Start();

			p_j_r = next->Prev()->Start()->Position();
			p_j_l = next->Twin()->Prev()->Start()->Position();

			double w = Cot(p_i->Position(), p_j_l, p_j->Position()) + Cot(p_i->Position(), p_j_r, p_j->Position());
			factor->AddElement(i, p_j->Index(), w*lamda);

			wj_sum += w;
			next = next->Prev()->Twin();
		} while (next && next != start);

		factor->SortMatrix();
		p_i->SetFlag(1);
		//cout << "factor:\n" << *factor << endl;
		factor->divideByRow(i, wj_sum);
	}
	// get result

	factor->Multiply(vx, vx1);
	factor->Multiply(vy, vy1);
	factor->Multiply(vz, vz1);
	
	/*cout <<"old:"<< vx <<"\n"<< vy << "\n" << vz << endl;
	cout << "factor:\n" << *factor << endl;
	cout <<"new:"<< vx1 << "\n" << vy1 << "\n" << vz1 << endl;*/

	// update coordinates
	for (size_t i = 0; i < v; i++) {
		vList[i]->SetPosition(Vector3d(vx1[i], vy1[i], vz1[i]));
	}

	//delete
	delete[]vx;
	delete[]vy;
	delete[]vz;
	delete[]vx1;
	delete[]vy1;
	delete[]vz1;
	delete factor;
}

void Mesh::ImplicitUmbrellaSmooth()
{
	/*************************/
	/* insert your code here */
	/*************************/
}


// helper for color ramp
Vector3d getColor(double curvature, double min, double max) {
	double r = (curvature-min)*0.7;
	double g = (max - curvature) * 0.01;
	double b = (curvature - min);
	Vector3d color = Vector3d(r, g, b);
	return color / color.L2Norm();
}
void Mesh::ComputeVertexCurvatures()
{
	/*************************/
	/* insert your code here */
	/*************************/
	double min = 100., max = 1.;
	VertexList l = vList;
	for (size_t i = 0; i < l.size(); i++) {
		// visit all vertices
		Vertex* p = l[i];
		int k = p->Valence();
		if (p->IsBoundary() || k < 3) {
			p->SetColor(Vector3d(0, 0, 0));
			p->SetFlag(1);
			continue;
		}

		double a = 0;
		OneRingVertex ring(p);
		Vector3d kn = Vector3d(0, 0, 0), p_j, p_j_l, p_j_r;
		HEdge* start = p->HalfEdge();
		HEdge* next = start;
		do{
			p_j = next->Next()->Start()->Position();
			p_j_r = next->Prev()->Start()->Position();
			p_j_l = next->Twin()->Prev()->Start()->Position();
			
			// update A
			a += Area(p->Position(), p_j, p_j_r);
			
			//update sum
			double cotSum = Cot(p->Position(), p_j_l, p_j) + Cot(p->Position(), p_j_r, p_j);
			kn += cotSum * (p_j - p->Position());
			
			// update next edge
			next = next->Prev()->Twin();
		} while (next && next != start);
		kn /= (-4 * a);
		double curv = kn.L2Norm();
		p->SetCurv(curv);
		//cout << "curv = " << curv << endl;
		if (curv > max) { max = curv; 
			//cout << kn << "\na=" << a << "\ncurv="<<curv<<endl;
		}
		if (curv < min) min = curv;		
	}
	//set color
	for (size_t i = 0; i < l.size(); i++) {
		// visit all vertices
		Vertex* p = l[i];
		if (p->Flag() ) {
			continue;
		}
		// set color
		p->SetColor(getColor(p->Curv(), min, max));
	}
	for (size_t i = 0; i < l.size(); i++) {
		l[i]->SetFlag(0);
	}
	//cout << "maxCurvature:"<<max << "\nminCurvature:"<< min<<endl;
	//0.4~394
}


