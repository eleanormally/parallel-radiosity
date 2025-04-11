#include <algorithm>
#include <cstring>
#include <iostream>

#include "argparser.h"
#include "face.h"
#include "kdtree.h"
#include "mesh.h"
#include "photon_mapping.h"
#include "primitive.h"
#include "raytracer.h"
#include "utils.h"

// ==========
// Clear/reset
void PhotonMapping::Clear() {
  // cleanup all the photons
  delete kdtree;
  kdtree = NULL;
}

// ========================================================================
// Recursively trace a single photon

void PhotonMapping::TracePhoton(const Vec3f& position, const Vec3f& direction,
                                const Vec3f& energy, int iter) {

  // ==============================================
  // ASSIGNMENT: IMPLEMENT RECURSIVE PHOTON TRACING
  // ==============================================
  Hit h;
  const Ray r(position, direction);
  if (!raytracer->CastRay(r, h, false)) {
    return;
  }
  Vec3f newPosition = position + (direction * (h.getT()));
  const Vec3f diffuseColor =
      h.getMaterial()->getDiffuseColor(h.get_s(), h.get_t());
  const Vec3f reflectiveColor = h.getMaterial()->getReflectiveColor();
  const Vec3f normal = h.getNormal();
  if (iter > 0) {
    Photon newPhoton(newPosition, direction, energy, iter);
    kdtree->AddPhoton(std::move(newPhoton));
  }

  Vec3f newDirection, newEnergy;
  if (reflectiveColor.Length() < EPSILON) {
    newDirection = RandomDiffuseDirection(normal);
    newEnergy = energy * diffuseColor;
  } else {
    newDirection = (direction - (2 * direction.Dot3(normal) * normal));
    newDirection.Normalize();
    newEnergy = energy * reflectiveColor;
  }
  if (iter >= GLOBAL_args->mesh_data->num_bounces) {
    return;
  }
  TracePhoton(newPosition, newDirection, newEnergy, iter + 1);
  // Trace the photon through the scene.  At each diffuse or
  // reflective bounce, store the photon in the kd tree.

  // One optimization is to *not* store the first bounce, since that
  // direct light can be efficiently computed using classic ray
  // tracing.
}

// ========================================================================
// Trace the specified number of photons through the scene

void PhotonMapping::TracePhotons() {

  // first, throw away any existing photons
  delete kdtree;

  // consruct a kdtree to store the photons
  BoundingBox* bb = mesh->getBoundingBox();
  Vec3f min = bb->getMin();
  Vec3f max = bb->getMax();
  Vec3f diff = max - min;
  min -= 0.001f * diff;
  max += 0.001f * diff;
  kdtree = new KDTree(BoundingBox(min, max));

  // photons emanate from the light sources
  const std::vector<Face*>& lights = mesh->getLights();

  // compute the total area of the lights
  float total_lights_area = 0;
  for (unsigned int i = 0; i < lights.size(); i++) {
    total_lights_area += lights[i]->getArea();
  }

  // shoot a constant number of photons per unit area of light source
  // (alternatively, this could be based on the total energy of each light)
  for (unsigned int i = 0; i < lights.size(); i++) {
    float my_area = lights[i]->getArea();
    int num =
        args->mesh_data->num_photons_to_shoot * my_area / total_lights_area;
    // the initial energy for this photon
    Vec3f energy =
        my_area / float(num) * lights[i]->getMaterial()->getEmittedColor();
    Vec3f normal = lights[i]->computeNormal();
    for (int j = 0; j < num; j++) {
      Vec3f start = lights[i]->RandomPoint();
      // the initial direction for this photon (for diffuse light sources)
      Vec3f direction = RandomDiffuseDirection(normal);
      TracePhoton(start, direction, energy, 0);
    }
  }
}

// ======================================================================

// helper function
bool closest_photon(const std::pair<Photon, float>& a,
                    const std::pair<Photon, float>& b) {
  return (a.second < b.second);
}

BoundingBox sizedBboxAboutPoint(const Vec3f& center, const float radius) {
  Vec3f delta(radius, radius, radius);
  return BoundingBox(center - delta, center + delta);
}

typedef std::pair<float, Vec3f> DistEnergy;
bool compareDistEnergy(const DistEnergy& a, const DistEnergy& b) {
  return a.first < b.first;
}
// ======================================================================
Vec3f PhotonMapping::GatherIndirect(const Vec3f& point, const Vec3f& /*normal*/,
                                    const Vec3f& /*direction_from*/) const {

  if (kdtree == NULL) {
    std::cout << "WARNING: Photons have not been traced throughout the scene."
              << std::endl;
    return Vec3f(0, 0, 0);
  }

  // ================================================================
  // ASSIGNMENT: GATHER THE INDIRECT ILLUMINATION FROM THE PHOTON MAP
  // ================================================================
  //binary searching bounding box
  size_t targetPhotonCount = GLOBAL_args->mesh_data->num_photons_to_collect;
#define MAX_RECURSES 8
  std::vector<Photon const*> photonsInRange;
  float currentRadius = 1;
  for (int _ = 0; _ < MAX_RECURSES; _++) {
    BoundingBox bb = sizedBboxAboutPoint(point, currentRadius / 2);
    photonsInRange.clear();
    kdtree->CollectPhotonsInBox(bb, photonsInRange);
    if (photonsInRange.size() >= targetPhotonCount) {
      break;
    }
    currentRadius *= 1.3;
  }
  if (photonsInRange.size() < targetPhotonCount) {
    return Vec3f();
  }

  std::vector<DistEnergy> squareDistEnergy(photonsInRange.size());
  for (size_t i = 0; i < photonsInRange.size(); i++) {
    const Vec3f& photonPos = photonsInRange[i]->getPosition();
    const Vec3f& photonEnergy = photonsInRange[i]->getEnergy();
    const float distToPoint = Vec3f::squareDist(photonPos, point);
    squareDistEnergy[i] = std::make_pair(distToPoint, photonEnergy);
  }
  std::sort(squareDistEnergy.begin(), squareDistEnergy.end(),
            compareDistEnergy);
  Vec3f totalEnergy;
  for (size_t i = 0; i < targetPhotonCount; i++) {
    totalEnergy += squareDistEnergy[i].second;
  }
  float boundingRadius = squareDistEnergy[targetPhotonCount - 1].first;
  totalEnergy /= boundingRadius * 3.14159;
  return totalEnergy;

  // collect the closest args->num_photons_to_collect photons
  // determine the radius that was necessary to collect that many photons
  // average the energy of those photons over that radius

  // return the color
}

// ======================================================================
// ======================================================================
// Helper functions to render the photons & kdtree

int PhotonMapping::triCount() const {
  int tri_count = 0;
  if (GLOBAL_args->mesh_data->render_kdtree == true && kdtree != NULL)
    tri_count += kdtree->numBoxes() * 12 * 12;
  if (GLOBAL_args->mesh_data->render_photon_directions == true &&
      kdtree != NULL)
    tri_count += kdtree->numPhotons() * 12;
  return tri_count;
}

int PhotonMapping::pointCount() const {
  if (GLOBAL_args->mesh_data->render_photons == false || kdtree == NULL)
    return 0;
  return kdtree->numPhotons();
}

// defined in raytree.cpp
void addBox(float*& current, Vec3f start, Vec3f end, Vec3f color, float width);

// ======================================================================

void packKDTree(const KDTree* kdtree, float*& current, int& count) {
  if (!kdtree->isLeaf()) {
    if (kdtree->getChild1() != NULL) {
      packKDTree(kdtree->getChild1(), current, count);
    }
    if (kdtree->getChild2() != NULL) {
      packKDTree(kdtree->getChild2(), current, count);
    }
  } else {

    Vec3f a = kdtree->getMin();
    Vec3f b = kdtree->getMax();

    Vec3f corners[8] = {Vec3f(a.x(), a.y(), a.z()), Vec3f(a.x(), a.y(), b.z()),
                        Vec3f(a.x(), b.y(), a.z()), Vec3f(a.x(), b.y(), b.z()),
                        Vec3f(b.x(), a.y(), a.z()), Vec3f(b.x(), a.y(), b.z()),
                        Vec3f(b.x(), b.y(), a.z()), Vec3f(b.x(), b.y(), b.z())};

    float width = 0.01 * (a - b).Length();

    addBox(current, corners[0], corners[1], Vec3f(1, 1, 0), width);
    addBox(current, corners[1], corners[3], Vec3f(1, 1, 0), width);
    addBox(current, corners[3], corners[2], Vec3f(1, 1, 0), width);
    addBox(current, corners[2], corners[0], Vec3f(1, 1, 0), width);

    addBox(current, corners[4], corners[5], Vec3f(1, 1, 0), width);
    addBox(current, corners[5], corners[7], Vec3f(1, 1, 0), width);
    addBox(current, corners[7], corners[6], Vec3f(1, 1, 0), width);
    addBox(current, corners[6], corners[4], Vec3f(1, 1, 0), width);

    addBox(current, corners[0], corners[4], Vec3f(1, 1, 0), width);
    addBox(current, corners[1], corners[5], Vec3f(1, 1, 0), width);
    addBox(current, corners[2], corners[6], Vec3f(1, 1, 0), width);
    addBox(current, corners[3], corners[7], Vec3f(1, 1, 0), width);

    count++;
  }
}

// ======================================================================

void packPhotons(const KDTree* kdtree, float*& current_points, int& count) {
  if (!kdtree->isLeaf()) {
    if (kdtree->getChild1() != NULL) {
      packPhotons(kdtree->getChild1(), current_points, count);
    }
    if (kdtree->getChild2() != NULL) {
      packPhotons(kdtree->getChild2(), current_points, count);
    }
  } else {
    for (unsigned int i = 0; i < kdtree->getPhotons().size(); i++) {
      const Photon& p = kdtree->getPhotons()[i];
      Vec3f v = p.getPosition();
      Vec3f color =
          p.getEnergy() * float(GLOBAL_args->mesh_data->num_photons_to_shoot);
      float12 t = {
          float(v.x()),     float(v.y()),     float(v.z()),     1, 0, 0, 0, 0,
          float(color.r()), float(color.g()), float(color.b()), 1};
      memcpy(current_points, &t, sizeof(float) * 12);
      current_points += 12;
      count++;
    }
  }
}

void packPhotonDirections(const KDTree* kdtree, float*& current, int& count) {
  if (!kdtree->isLeaf()) {
    if (kdtree->getChild1() != NULL) {
      packPhotonDirections(kdtree->getChild1(), current, count);
    }
    if (kdtree->getChild2() != NULL) {
      packPhotonDirections(kdtree->getChild2(), current, count);
    }
  } else {
    for (unsigned int i = 0; i < kdtree->getPhotons().size(); i++) {
      const Photon& p = kdtree->getPhotons()[i];
      Vec3f v = p.getPosition();
      Vec3f v2 = p.getPosition() - p.getDirectionFrom() * 0.5;
      Vec3f color =
          p.getEnergy() * float(GLOBAL_args->mesh_data->num_photons_to_shoot);
      float width = 0.01;
      addBox(current, v, v2, color, width);
      count++;
    }
  }
}

// ======================================================================

void PhotonMapping::packMesh(float*& current, float*& current_points) {

  // the photons
  if (GLOBAL_args->mesh_data->render_photons && kdtree != NULL) {
    int count = 0;
    packPhotons(kdtree, current_points, count);
    assert(count == kdtree->numPhotons());
  }
  // photon directions
  if (GLOBAL_args->mesh_data->render_photon_directions && kdtree != NULL) {
    int count = 0;
    packPhotonDirections(kdtree, current, count);
    assert(count == kdtree->numPhotons());
  }

  // the wireframe kdtree
  if (GLOBAL_args->mesh_data->render_kdtree && kdtree != NULL) {
    int count = 0;
    packKDTree(kdtree, current, count);
    assert(count == kdtree->numBoxes());
  }
}

// ======================================================================
