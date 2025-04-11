#include "radiosity.h"
#include "face.h"
#include "mesh.h"
#include "raytracer.h"
#include "raytree.h"
#include "sphere.h"
#include "utils.h"
#include "vectors.h"

// ================================================================
// CONSTRUCTOR & DESTRUCTOR
// ================================================================
Radiosity::Radiosity(Mesh* m, ArgParser* a) {
  mesh = m;
  args = a;
  num_faces = -1;
  formfactors = NULL;
  area = NULL;
  undistributed = NULL;
  absorbed = NULL;
  radiance = NULL;
  max_undistributed_patch = -1;
  total_area = -1;
  Reset();
}

Radiosity::~Radiosity() {
  Cleanup();
}

void Radiosity::Cleanup() {
  delete[] formfactors;
  delete[] area;
  delete[] undistributed;
  delete[] absorbed;
  delete[] radiance;
  num_faces = -1;
  formfactors = NULL;
  area = NULL;
  undistributed = NULL;
  absorbed = NULL;
  radiance = NULL;
  max_undistributed_patch = -1;
  total_area = -1;
}

void Radiosity::Reset() {
  delete[] area;
  delete[] undistributed;
  delete[] absorbed;
  delete[] radiance;

  // create and fill the data structures
  num_faces = mesh->numFaces();
  area = new float[num_faces];
  undistributed = new Vec3f[num_faces];
  absorbed = new Vec3f[num_faces];
  radiance = new Vec3f[num_faces];
  for (int i = 0; i < num_faces; i++) {
    Face* f = mesh->getFace(i);
    f->setRadiosityPatchIndex(i);
    setArea(i, f->getArea());
    Vec3f emit = f->getMaterial()->getEmittedColor();
    setUndistributed(i, emit);
    setAbsorbed(i, Vec3f(0, 0, 0));
    setRadiance(i, emit);
  }

  // find the patch with the most undistributed energy
  findMaxUndistributed();
}

// =======================================================================================
// =======================================================================================

void Radiosity::findMaxUndistributed() {
  // find the patch with the most undistributed energy
  // don't forget that the patches may have different sizes!
  max_undistributed_patch = -1;
  total_undistributed = 0;
  total_area = 0;
  float max = -1;
  for (int i = 0; i < num_faces; i++) {
    float m = getUndistributed(i).Length() * getArea(i);
    total_undistributed += m;
    total_area += getArea(i);
    if (max < m) {
      max = m;
      max_undistributed_patch = i;
    }
  }
  assert(max_undistributed_patch >= 0 && max_undistributed_patch < num_faces);
}

void Radiosity::ComputeFormFactors() {
  /*assert(formfactors == NULL);*/
  if (formfactors != nullptr) {
    delete[] formfactors;
  }
  assert(num_faces > 0);
  const int samples = GLOBAL_args->mesh_data->num_form_factor_samples;
  const int shadowSamples = GLOBAL_args->mesh_data->num_shadow_samples;
  formfactors = new float[num_faces * num_faces];
  //calculating form factors with shadows
  for (int i = 0; i < num_faces; i++) {
    for (int j = 0; j < num_faces; j++) {
      int         idx = (i * num_faces) + j;
      const Face* iFace = mesh->getFace(i);
      const Face* jFace = mesh->getFace(j);
      if (j < i && false) {
        formfactors[idx] = formfactors[(j * num_faces) + i] / jFace->getArea() *
                           iFace->getArea();
        continue;
      }
      if (i == j) {
        formfactors[idx] = 0;
        continue;
      }
      const Vec3f iNorm = iFace->computeNormal();
      const Vec3f jNorm = jFace->computeNormal();
      double      totalFlux = 0.0;
      for (int _ = 0; _ < samples; _++) {
        const Vec3f  iPoint = iFace->RandomPoint();
        const Vec3f  jPoint = jFace->RandomPoint();
        const Vec3f  iToJ = jPoint - iPoint;
        const Vec3f  jToI = iPoint - jPoint;
        const double iCos = iNorm.Dot3(iToJ) / iToJ.Length();
        const double jCos = jNorm.Dot3(jToI) / jToI.Length();
        totalFlux += fabs(iCos * jCos) / (jToI.Length() * jToI.Length());
      }
      float totalVisible = 0.0;
      for (int _ = 0; _ < shadowSamples; _++) {
        Vec3f iPoint, jPoint;
        if (shadowSamples == 1) {
          iPoint = iFace->computeCentroid();
          jPoint = jFace->computeCentroid();
        } else {
          iPoint = iFace->RandomPoint();
          jPoint = jFace->RandomPoint();
        }
        Vec3f        iToJ = jPoint - iPoint;
        const double iToJLength = iToJ.Length();
        iToJ.Normalize();
        const Ray r(iPoint, iToJ);
        Hit       h;
        raytracer->CastRay(r, h, true);
        if (fabs(h.getT() - iToJLength) < EPSILON) {
          totalVisible += 1.0;
        }
      }
      if (shadowSamples > 0) {
        totalVisible /= shadowSamples;
        totalFlux *= totalVisible;
      }
      totalFlux /= samples;
      totalFlux /= 3.14159;
      totalFlux /= iFace->getArea();
      formfactors[idx] = totalFlux;
    }
  }

  //normalizing each row so all light hits something
  for (int i = 0; i < num_faces; i++) {
    double rowFlux = 0;
    for (int j = 0; j < num_faces; j++) {
      int idx = (i * num_faces) + j;
      rowFlux += formfactors[idx];
    }
    if (rowFlux == 0)
      continue;
    for (int j = 0; j < num_faces; j++) {
      int idx = (i * num_faces) + j;
      formfactors[idx] /= rowFlux;
    }
  }

  // =====================================
  // ASSIGNMENT:  COMPUTE THE FORM FACTORS
  // =====================================
}

// ================================================================
// ================================================================

void Radiosity::distributeExtraRadiosityToPatch(int source, int receiver) {
  const float proportionToReceiver = getFormFactor(source, receiver);
  const Vec3f undistributedSource = getUndistributed(source);
  const Vec3f incomingRadiosity = undistributedSource * proportionToReceiver;

  const Material* receiverMaterial = mesh->getFace(receiver)->getMaterial();
  const Vec3f     proportionReflected = receiverMaterial->getDiffuseColor();
  const Vec3f     proportionAbsorbed = Vec3f(1, 1, 1) - proportionReflected;

  absorbed[receiver] += proportionAbsorbed * incomingRadiosity;
  undistributed[receiver] += proportionReflected * incomingRadiosity;
  radiance[receiver] += proportionReflected * incomingRadiosity;
}

float Radiosity::Iterate() {
  if (formfactors == NULL) {
    ComputeFormFactors();
    for (int i = 0; i < num_faces; i++) {
      for (int j = 0; j < num_faces; j++) {
        std::cout << getFormFactor(i, j) << ", ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "diffuse:\n";
    for (int i = 0; i < num_faces; i++) {
      std::cout << mesh->getFace(i)->getMaterial()->getDiffuseColor() << ", ";
    }
    std::cout << "\n\nemmitence:\n";
    for (int i = 0; i < num_faces; i++) {
      std::cout << mesh->getFace(i)->getMaterial()->getEmittedColor() << ", ";
    }
    std::cout << std::endl;
  }
  assert(formfactors != NULL);

  // ==========================================
  // ASSIGNMENT:  IMPLEMENT RADIOSITY ALGORITHM
  // ==========================================

  findMaxUndistributed();

  int source = max_undistributed_patch;
  for (int receiving = 0; receiving < num_faces; receiving++) {
    distributeExtraRadiosityToPatch(source, receiving);
  }
  undistributed[source] = Vec3f(0, 0, 0);

  float totalUndistributed = 0;
  for (int i = 0; i < num_faces; i++) {
    totalUndistributed += undistributed[i].Length();
  }
  return totalUndistributed;
}

// =======================================================================================
// HELPER FUNCTIONS FOR RENDERING
// =======================================================================================

// for interpolation
void CollectFacesWithVertex(Vertex* have, Face* f, std::vector<Face*>& faces) {
  for (unsigned int i = 0; i < faces.size(); i++) {
    if (faces[i] == f)
      return;
  }
  if (have != (*f)[0] && have != (*f)[1] && have != (*f)[2] && have != (*f)[3])
    return;
  faces.push_back(f);
  for (int i = 0; i < 4; i++) {
    Edge* ea = f->getEdge()->getOpposite();
    Edge* eb = f->getEdge()->getNext()->getOpposite();
    Edge* ec = f->getEdge()->getNext()->getNext()->getOpposite();
    Edge* ed = f->getEdge()->getNext()->getNext()->getNext()->getOpposite();
    if (ea != NULL)
      CollectFacesWithVertex(have, ea->getFace(), faces);
    if (eb != NULL)
      CollectFacesWithVertex(have, eb->getFace(), faces);
    if (ec != NULL)
      CollectFacesWithVertex(have, ec->getFace(), faces);
    if (ed != NULL)
      CollectFacesWithVertex(have, ed->getFace(), faces);
  }
}

// different visualization modes
Vec3f Radiosity::setupHelperForColor(Face* f, int i, int j) {
  assert(mesh->getFace(i) == f);
  assert(j >= 0 && j < 4);
  if (args->mesh_data->render_mode == RENDER_MATERIALS) {
    return f->getMaterial()->getDiffuseColor();
  } else if (args->mesh_data->render_mode == RENDER_RADIANCE &&
             args->mesh_data->interpolate == true) {
    std::vector<Face*> faces;
    CollectFacesWithVertex((*f)[j], f, faces);
    float total = 0;
    Vec3f color = Vec3f(0, 0, 0);
    Vec3f normal = f->computeNormal();
    for (unsigned int i = 0; i < faces.size(); i++) {
      Vec3f normal2 = faces[i]->computeNormal();
      float area = faces[i]->getArea();
      if (normal.Dot3(normal2) < 0.5)
        continue;
      assert(area > 0);
      total += area;
      color += float(area) * getRadiance(faces[i]->getRadiosityPatchIndex());
    }
    assert(total > 0);
    color /= total;
    return color;
  } else if (args->mesh_data->render_mode == RENDER_LIGHTS) {
    return f->getMaterial()->getEmittedColor();
  } else if (args->mesh_data->render_mode == RENDER_UNDISTRIBUTED) {
    return getUndistributed(i);
  } else if (args->mesh_data->render_mode == RENDER_ABSORBED) {
    return getAbsorbed(i);
  } else if (args->mesh_data->render_mode == RENDER_RADIANCE) {
    return getRadiance(i);
  } else if (args->mesh_data->render_mode == RENDER_FORM_FACTORS) {
    if (formfactors == NULL)
      ComputeFormFactors();
    float scale = 0.2 * total_area / getArea(i);
    float factor = scale * getFormFactor(max_undistributed_patch, i);
    return Vec3f(factor, factor, factor);
  } else {
    assert(0);
  }
  exit(0);
}

// =======================================================================================

int Radiosity::triCount() {
  return 12 * num_faces;
}

void Radiosity::packMesh(float*& current) {

  for (int i = 0; i < num_faces; i++) {
    Face* f = mesh->getFace(i);
    Vec3f normal = f->computeNormal();

    //double avg_s = 0;
    //double avg_t = 0;

    // wireframe is normally black, except when it's the special
    // patch, then the wireframe is red
    Vec3f wireframe_color(0, 0, 0);
    if (args->mesh_data->render_mode == RENDER_FORM_FACTORS &&
        i == max_undistributed_patch) {
      wireframe_color = Vec3f(1, 0, 0);
    }

    // 4 corner vertices
    Vec3f a_pos = ((*f)[0])->get();
    Vec3f a_color = setupHelperForColor(f, i, 0);
    a_color = Vec3f(linear_to_srgb(a_color.r()), linear_to_srgb(a_color.g()),
                    linear_to_srgb(a_color.b()));
    Vec3f b_pos = ((*f)[1])->get();
    Vec3f b_color = setupHelperForColor(f, i, 1);
    b_color = Vec3f(linear_to_srgb(b_color.r()), linear_to_srgb(b_color.g()),
                    linear_to_srgb(b_color.b()));
    Vec3f c_pos = ((*f)[2])->get();
    Vec3f c_color = setupHelperForColor(f, i, 2);
    c_color = Vec3f(linear_to_srgb(c_color.r()), linear_to_srgb(c_color.g()),
                    linear_to_srgb(c_color.b()));
    Vec3f d_pos = ((*f)[3])->get();
    Vec3f d_color = setupHelperForColor(f, i, 3);
    d_color = Vec3f(linear_to_srgb(d_color.r()), linear_to_srgb(d_color.g()),
                    linear_to_srgb(d_color.b()));

    Vec3f avg_color = 0.25f * (a_color + b_color + c_color + d_color);

    // the centroid (for wireframe rendering)
    Vec3f centroid = f->computeCentroid();

    AddWireFrameTriangle(current, a_pos, b_pos, centroid, normal, normal,
                         normal, wireframe_color, a_color, b_color, avg_color);
    AddWireFrameTriangle(current, b_pos, c_pos, centroid, normal, normal,
                         normal, wireframe_color, b_color, c_color, avg_color);
    AddWireFrameTriangle(current, c_pos, d_pos, centroid, normal, normal,
                         normal, wireframe_color, c_color, d_color, avg_color);
    AddWireFrameTriangle(current, d_pos, a_pos, centroid, normal, normal,
                         normal, wireframe_color, d_color, a_color, avg_color);
  }
}

// =======================================================================================
