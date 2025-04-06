#include "sphere.h"
#include "argparser.h"
#include "hit.h"
#include "material.h"
#include "mesh.h"
#include "ray.h"
#include "utils.h"
#include "vertex.h"

// ====================================================================
// ====================================================================

bool Sphere::intersect(const Ray& r, Hit& h) const {

  //representing the ray as an offset of the sphere so the sphere
  //pretends it is at the origin

  Vec3f origin = r.getOrigin() - center;
  const double a = 1;
  const double b = 2 * r.getDirection().Dot3(origin);
  const double c = origin.Dot3(origin) - (radius * radius);

  const double discriminantSquared = (b * b) - (4 * a * c);
  if (discriminantSquared <= 0) {
    return false;
  }
  const double discriminant = sqrt(discriminantSquared);

  double tSolution1 = (-1 * b) + discriminant;
  double tSolution2 = (-1 * b) - discriminant;
  double tSolution;
  if (tSolution2 < 0) {
    tSolution = tSolution1;
  } else {
    tSolution = std::min(tSolution1, tSolution2);
  }
  tSolution /= 2 * a;
  if (tSolution < EPSILON)
    return false;
  if (h.getT() < tSolution) {
    return true;
  }
  Vec3f intersection = r.pointAtParameter(tSolution);
  Vec3f normal = intersection - center;
  normal.Normalize();
  h.set(tSolution, material, normal);
  // return true if the sphere was intersected, and update the hit
  // data structure to contain the value of t for the ray at the
  // intersection point, the material, and the normal

  return true;
}

// ====================================================================
// ====================================================================

// helper function to place a grid of points on the sphere
Vec3f ComputeSpherePoint(float s, float t, const Vec3f center, float radius) {
  float angle = 2 * M_PI * s;
  float y = -cos(M_PI * t);
  float factor = sqrt(1 - y * y);
  float x = factor * cos(angle);
  float z = factor * -sin(angle);
  Vec3f answer = Vec3f(x, y, z);
  answer *= radius;
  answer += center;
  return answer;
}

void Sphere::addRasterizedFaces(Mesh* m, ArgParser* args) {

  // and convert it into quad patches for radiosity
  int h = args->mesh_data->sphere_horiz;
  int v = args->mesh_data->sphere_vert;
  assert(h % 2 == 0);
  int i, j;
  int va, vb, vc, vd;
  Vertex *a, *b, *c, *d;
  int offset = m->numVertices();  //vertices.size();

  // place vertices
  m->addVertex(center + radius * Vec3f(0, -1, 0));  // bottom
  for (j = 1; j < v; j++) {                         // middle
    for (i = 0; i < h; i++) {
      float s = i / float(h);
      float t = j / float(v);
      m->addVertex(ComputeSpherePoint(s, t, center, radius));
    }
  }
  m->addVertex(center + radius * Vec3f(0, 1, 0));  // top

  // the middle patches
  for (j = 1; j < v - 1; j++) {
    for (i = 0; i < h; i++) {
      va = 1 + i + h * (j - 1);
      vb = 1 + (i + 1) % h + h * (j - 1);
      vc = 1 + i + h * (j);
      vd = 1 + (i + 1) % h + h * (j);
      a = m->getVertex(offset + va);
      b = m->getVertex(offset + vb);
      c = m->getVertex(offset + vc);
      d = m->getVertex(offset + vd);
      m->addRasterizedPrimitiveFace(a, b, d, c, material);
    }
  }

  for (i = 0; i < h; i += 2) {
    // the bottom patches
    va = 0;
    vb = 1 + i;
    vc = 1 + (i + 1) % h;
    vd = 1 + (i + 2) % h;
    a = m->getVertex(offset + va);
    b = m->getVertex(offset + vb);
    c = m->getVertex(offset + vc);
    d = m->getVertex(offset + vd);
    m->addRasterizedPrimitiveFace(d, c, b, a, material);
    // the top patches
    va = 1 + h * (v - 1);
    vb = 1 + i + h * (v - 2);
    vc = 1 + (i + 1) % h + h * (v - 2);
    vd = 1 + (i + 2) % h + h * (v - 2);
    a = m->getVertex(offset + va);
    b = m->getVertex(offset + vb);
    c = m->getVertex(offset + vc);
    d = m->getVertex(offset + vd);
    m->addRasterizedPrimitiveFace(b, c, d, a, material);
  }
}
