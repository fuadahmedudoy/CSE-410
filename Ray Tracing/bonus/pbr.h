#pragma once

#include <bits/stdc++.h>
#include <GL/glut.h>
#include "bitmap_image.hpp"

#define EPSILON 1e-7
#define PI 3.14159265358979323846

extern bitmap_image image;
extern int recursionLevel;

double whiteColor[] = {1, 1, 1};
double blackColor[] = {0, 0, 0};

// Color class for texture sampling
class Color {
public:
    double r, g, b;
    
    Color() : r(0), g(0), b(0) {}
    
    Color(double r, double g, double b) : r(r), g(g), b(b) {}
    
    double* toArray() {
        static double arr[3];
        arr[0] = r;
        arr[1] = g;
        arr[2] = b;
        return arr;
    }
    
    Color operator*(double scalar) const {
        return Color(r * scalar, g * scalar, b * scalar);
    }
    
    Color operator+(const Color& other) const {
        return Color(r + other.r, g + other.g, b + other.b);
    }
    
    // Add Color * Color operator
    Color operator*(const Color& other) const {
        return Color(r * other.r, g * other.g, b * other.b);
    }
};

// PBR Material Properties
struct PBRMaterial {
    Color albedo;        // Base color (diffuse reflectance)
    double metallic;     // Metallic factor (0 = dielectric, 1 = metal)
    double roughness;    // Surface roughness (0 = smooth, 1 = rough)
    double ao;           // Ambient occlusion
    double reflectance;  // Reflectance at normal incidence (F0)
    
    PBRMaterial() : albedo(0.5, 0.5, 0.5), metallic(0.0), roughness(0.5), ao(1.0), reflectance(0.04) {}
    
    PBRMaterial(const Color& albedo, double metallic, double roughness, double ao) 
        : albedo(albedo), metallic(metallic), roughness(roughness), ao(ao), reflectance(0.04) {}
};

// Global texture variables
extern unsigned char* textureData;
extern int textureWidth;
extern int textureHeight;
extern int textureChannels;

using namespace std;

class Vector3D
{
public:
    double x, y, z;

    Vector3D()
    {
        x = 0;
        y = 0;
        z = 0;
    }

    Vector3D(const Vector3D &v)
    {
        x = v.x;
        y = v.y;
        z = v.z;
    }

    Vector3D(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    Vector3D &operator=(const Vector3D &v)
    {
        if (this != &v)
        {
            x = v.x;
            y = v.y;
            z = v.z;
        }
        return *this;
    }

    Vector3D operator+(const Vector3D &v)
    {
        return Vector3D(x + v.x, y + v.y, z + v.z);
    }

    Vector3D operator-(const Vector3D &v)
    {
        return Vector3D(x - v.x, y - v.y, z - v.z);
    }

    Vector3D operator+(double s)
    {
        return Vector3D(x + s, y + s, z + s);
    }

    Vector3D operator*(double s)
    {
        return Vector3D(x * s, y * s, z * s);
    }

    Vector3D operator/(double s)
    {
        assert(fabs(s) > EPSILON);
        return Vector3D(x / s, y / s, z / s);
    }

    // Change dot to be const method
    double dot(const Vector3D &v) const
    {
        return x * v.x + y * v.y + z * v.z;
    }

    Vector3D cross(const Vector3D &v) const
    {
        return Vector3D(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x);
    }

    double length() const
    {
        return sqrt(x * x + y * y + z * z);
    }

    Vector3D &normalize()
    {
        double l = length();
        assert(fabs(l) > EPSILON);
        x /= l;
        y /= l;
        z /= l;
        return *this;
    }

    // Add a const version that returns a new normalized vector without modifying this one
    Vector3D normalized() const
    {
        double l = length();
        assert(fabs(l) > EPSILON);
        return Vector3D(x / l, y / l, z / l);
    }

    Vector3D rotateAroundAxis(double angle, Vector3D &axis) const
    {
        double s = sin(angle);
        double c = cos(angle);
        double x = this->x;
        double y = this->y;
        double z = this->z;
        double u = axis.x;
        double v = axis.y;
        double w = axis.z;
        double xPrime = (u * (u * x + v * y + w * z) * (1 - c) + x * c + (-w * y + v * z) * s);
        double yPrime = (v * (u * x + v * y + w * z) * (1 - c) + y * c + (w * x - u * z) * s);
        double zPrime = (w * (u * x + v * y + w * z) * (1 - c) + z * c + (-v * x + u * y) * s);
        return Vector3D(xPrime, yPrime, zPrime);
    }
};

class Ray
{
public:
    Vector3D start, dir;
    Ray(const Vector3D &start, const Vector3D &dir)
    {
        this->start = start;
        this->dir = dir.normalized(); // Use normalized() instead of normalize()
    }
};

class PointLight
{
public:
    Vector3D position;
    double color[3];
    double intensity;  // Added for PBR

    PointLight(const Vector3D &position, const double *color, double intensity = 1.0)
    {
        this->position = position;
        copy(color, color + 3, this->color);
        this->intensity = intensity;
    }
};

class SpotLight
{
public:
    Vector3D position;
    Vector3D direction;
    double cutOffAngle;
    double color[3];
    double intensity;  // Added for PBR
    
    SpotLight(const Vector3D &position, const Vector3D &direction, double cutOffAngle, double *color, double intensity = 1.0)
    {
        this->position = position;
        this->direction = direction;
        this->direction.normalize();
        this->cutOffAngle = cutOffAngle;
        copy(color, color + 3, this->color);
        this->intensity = intensity;
    }
};

class Object;
extern vector<Object *> objects;
extern vector<PointLight> pointLights;
extern vector<SpotLight> spotLights;

// PBR helper functions
inline double DistributionGGX(const Vector3D& N, const Vector3D& H, double roughness) {
    double a = roughness * roughness;
    double a2 = a * a;
    double NdotH = max(0.0, N.dot(H));
    double NdotH2 = NdotH * NdotH;
    
    double nom = a2;
    double denom = (NdotH2 * (a2 - 1.0) + 1.0);
    denom = PI * denom * denom;
    
    return nom / max(denom, EPSILON);
}

inline double GeometrySchlickGGX(double NdotV, double roughness) {
    double r = (roughness + 1.0);
    double k = (r * r) / 8.0;
    
    double nom = NdotV;
    double denom = NdotV * (1.0 - k) + k;
    
    return nom / max(denom, EPSILON);
}

inline double GeometrySmith(const Vector3D& N, const Vector3D& V, const Vector3D& L, double roughness) {
    double NdotV = max(0.0, N.dot(V));
    double NdotL = max(0.0, N.dot(L));
    double ggx2 = GeometrySchlickGGX(NdotV, roughness);
    double ggx1 = GeometrySchlickGGX(NdotL, roughness);
    
    return ggx1 * ggx2;
}

inline Color fresnelSchlick(double cosTheta, const Color& F0) {
    return Color(
        F0.r + (1.0 - F0.r) * pow(max(1.0 - cosTheta, 0.0), 5.0),
        F0.g + (1.0 - F0.g) * pow(max(1.0 - cosTheta, 0.0), 5.0),
        F0.b + (1.0 - F0.b) * pow(max(1.0 - cosTheta, 0.0), 5.0)
    );
}

class Object
{
public:
    Vector3D referencePoint;
    int shine;  // Legacy parameter
    double height, width, length;
    double color[3];
    double coEfficients[4];  // Legacy parameters
    PBRMaterial pbrMaterial;  // New PBR material properties

    virtual void setShine(int shine)
    {
        this->shine = shine;
    };

    virtual void setColor(double *color)
    {
        copy(color, color + 3, this->color);
        // Also set the albedo for PBR
        pbrMaterial.albedo = Color(color[0], color[1], color[2]);
    };

    virtual void setCoEfficients(double *coEfficients)
    {
        copy(coEfficients, coEfficients + 4, this->coEfficients);
        
        // Map legacy coefficients to PBR parameters
        // coEfficients[0] = ambient, [1] = diffuse, [2] = specular, [3] = reflection
        pbrMaterial.ao = coEfficients[0] * 2.0;  // Ambient to AO
        pbrMaterial.roughness = 1.0 - sqrt(coEfficients[2]);  // Inverse of specular
        pbrMaterial.metallic = coEfficients[3];  // Reflection to metallic
    };
    
    // New PBR material setters
    virtual void setPBRMaterial(const PBRMaterial& material) {
        this->pbrMaterial = material;
    }
    
    virtual void setMetallic(double metallic) {
        pbrMaterial.metallic = max(0.0, min(1.0, metallic));
    }
    
    virtual void setRoughness(double roughness) {
        pbrMaterial.roughness = max(0.0, min(1.0, roughness));
    }

    virtual double *getColor(Vector3D &point)
    {
        return this->color;
    }

    virtual Ray getReflectedRay(const Ray &incidentRay, const Vector3D &intersectionPoint)
    {
        Vector3D normal = this->getNormal(incidentRay, intersectionPoint);
        Vector3D reflectedDir = (normal * (2 * normal.dot(incidentRay.dir))) * -1.0 + incidentRay.dir;
        Vector3D reflectedStart = intersectionPoint;
        reflectedStart = reflectedStart + reflectedDir * EPSILON * 2.0;
        return Ray(reflectedStart, reflectedDir); // Ray constructor will normalize the direction
    }

    virtual double intersect(Ray &ray, double *color, int level)
    {
        double tMin = this->intersect(ray);
        if (level == 0)
        {
            return tMin;
        }

        Vector3D intersectionPoint = ray.start + ray.dir * tMin;
        double *intersectionPointColor = this->getColor(intersectionPoint);
        
        // Initialize color with ambient term
        color[0] = intersectionPointColor[0] * this->coEfficients[0] * pbrMaterial.ao;
        color[1] = intersectionPointColor[1] * this->coEfficients[0] * pbrMaterial.ao;
        color[2] = intersectionPointColor[2] * this->coEfficients[0] * pbrMaterial.ao;

        // PBR parameters
        Vector3D N = this->getNormal(ray, intersectionPoint);
        Vector3D V = (ray.start - intersectionPoint).normalized(); // Use normalized() instead of normalize()
        
        // Calculate F0 (surface reflection at zero incidence)
        Color F0(pbrMaterial.reflectance, pbrMaterial.reflectance, pbrMaterial.reflectance);
        
        // If material is metallic, F0 = albedo
        if (pbrMaterial.metallic > 0.0) {
            F0 = Color(
                F0.r * (1.0 - pbrMaterial.metallic) + pbrMaterial.albedo.r * pbrMaterial.metallic,
                F0.g * (1.0 - pbrMaterial.metallic) + pbrMaterial.albedo.g * pbrMaterial.metallic,
                F0.b * (1.0 - pbrMaterial.metallic) + pbrMaterial.albedo.b * pbrMaterial.metallic
            );
        }

        // Process point lights
        for (PointLight pointLight : pointLights)
        {
            Vector3D L = (pointLight.position - intersectionPoint).normalized(); // Use normalized()
            Vector3D H = (V + L).normalized(); // Use normalized()
            
            // Check if point is in shadow
            Vector3D originL = pointLight.position;
            Vector3D directionL = (intersectionPoint - originL).normalized(); // Use normalized()
            Ray lightRay(originL, directionL);
            bool isObscured = false;
            double len = (intersectionPoint - originL).length();

            if (len < EPSILON) {
                continue;
            }

            for (Object *object : objects)
            {
                double t = object->intersect(lightRay);
                if (t > EPSILON && (t - len) < -EPSILON)
                {
                    isObscured = true;
                    break;
                }
            }
            
            if (!isObscured)
            {
                // Calculate PBR lighting
                double distance = len;
                double attenuation = 1.0 / (1.0 + 0.09 * distance + 0.032 * distance * distance);
                double radiance = pointLight.intensity * attenuation;
                
                // Cook-Torrance BRDF
                double NdotV = max(0.0, N.dot(V));
                double NdotL = max(0.0, N.dot(L));
                
                double NDF = DistributionGGX(N, H, pbrMaterial.roughness);
                double G = GeometrySmith(N, V, L, pbrMaterial.roughness);
                Color F = fresnelSchlick(max(0.0, H.dot(V)), F0);
                
                // Calculate specular component
                Color numerator(
                    NDF * G * F.r,
                    NDF * G * F.g,
                    NDF * G * F.b
                );
                double denominator = 4.0 * NdotV * NdotL + EPSILON;
                Color specular(
                    numerator.r / denominator,
                    numerator.g / denominator,
                    numerator.b / denominator
                );
                
                // Calculate diffuse component (energy conservation)
                Color kS(F.r, F.g, F.b);
                Color kD(
                    (1.0 - kS.r) * (1.0 - pbrMaterial.metallic),
                    (1.0 - kS.g) * (1.0 - pbrMaterial.metallic),
                    (1.0 - kS.b) * (1.0 - pbrMaterial.metallic)
                );
                
                // Combine diffuse and specular
                Color Lo = (kD * Color(
                    pbrMaterial.albedo.r / PI,
                    pbrMaterial.albedo.g / PI,
                    pbrMaterial.albedo.b / PI
                ) + specular) * radiance * NdotL;
                
                // Add to final color
                color[0] += Lo.r * pointLight.color[0];
                color[1] += Lo.g * pointLight.color[1];
                color[2] += Lo.b * pointLight.color[2];
            }
        }

        // Process spot lights with PBR
        for (SpotLight spotLight : spotLights)
        {
            Vector3D originL = spotLight.position;
            Vector3D directionL = (intersectionPoint - originL).normalized(); // Use normalized()
            double cosBeta = directionL.dot(spotLight.direction);
            double beta = acos(cosBeta) * 180.0 / M_PI;
            double Xm = 0;
            double len = (intersectionPoint - originL).length();
            
            if (beta > spotLight.cutOffAngle || len < EPSILON)
            {
                continue;
            }
            else
            {
                double X = cos(spotLight.cutOffAngle * M_PI / 180.0);
                assert(X > 0);
                double d = 1.0 / (1.0 - X);
                if (X <= cosBeta && cosBeta <= 1.0)
                {
                    Xm = 1.0 - (1.0 - cosBeta) * d;
                }
            }
            
            Ray lightRay(originL, directionL);
            bool isObscured = false;
            
            for (Object *object : objects)
            {
                double t = object->intersect(lightRay);
                if (t > 0 && (t - len) < -EPSILON)
                {
                    isObscured = true;
                    break;
                }
            }
            
            if (!isObscured)
            {
                // PBR calculations for spot light
                Vector3D L = (spotLight.position - intersectionPoint).normalized(); // Use normalized()
                Vector3D H = (V + L).normalized(); // Use normalized()
                
                double distance = len;
                double attenuation = 1.0 / (1.0 + 0.09 * distance + 0.032 * distance * distance);
                double radiance = spotLight.intensity * attenuation * Xm;
                
                // Cook-Torrance BRDF
                double NdotV = max(0.0, N.dot(V));
                double NdotL = max(0.0, N.dot(L));
                
                double NDF = DistributionGGX(N, H, pbrMaterial.roughness);
                double G = GeometrySmith(N, V, L, pbrMaterial.roughness);
                Color F = fresnelSchlick(max(0.0, H.dot(V)), F0);
                
                // Calculate specular component
                Color numerator(
                    NDF * G * F.r,
                    NDF * G * F.g,
                    NDF * G * F.b
                );
                double denominator = 4.0 * NdotV * NdotL + EPSILON;
                Color specular(
                    numerator.r / denominator,
                    numerator.g / denominator,
                    numerator.b / denominator
                );
                
                // Calculate diffuse component (energy conservation)
                Color kS(F.r, F.g, F.b);
                Color kD(
                    (1.0 - kS.r) * (1.0 - pbrMaterial.metallic),
                    (1.0 - kS.g) * (1.0 - pbrMaterial.metallic),
                    (1.0 - kS.b) * (1.0 - pbrMaterial.metallic)
                );
                
                // Combine diffuse and specular
                Color Lo = (kD * Color(
                    pbrMaterial.albedo.r / PI,
                    pbrMaterial.albedo.g / PI,
                    pbrMaterial.albedo.b / PI
                ) + specular) * radiance * NdotL;
                
                // Add to final color
                color[0] += Lo.r * spotLight.color[0];
                color[1] += Lo.g * spotLight.color[1];
                color[2] += Lo.b * spotLight.color[2];
            }
        }

        // Handle reflections
        if (level >= recursionLevel)
        {
            return tMin;
        }

        Ray reflectedRay = this->getReflectedRay(ray, intersectionPoint);
        tMin = DBL_MAX;
        Object *nearestObject = nullptr;

        for (Object *object : objects)
        {
            double t = object->intersect(reflectedRay);
            if (t > 0 && t < tMin)
            {
                tMin = t;
                nearestObject = object;
            }
        }

        if (nearestObject != nullptr)
        {
            double *colorReflected = new double[3];
            tMin = nearestObject->intersect(reflectedRay, colorReflected, level + 1);
            
            // Use Fresnel factor for reflection strength
            double reflStrength = pbrMaterial.metallic + (1.0 - pbrMaterial.metallic) * 
                                 pow(1.0 - max(0.0, N.dot(V)), 5.0);
            
            color[0] += colorReflected[0] * reflStrength;
            color[1] += colorReflected[1] * reflStrength;
            color[2] += colorReflected[2] * reflStrength;
            delete[] colorReflected;
        }

        return tMin;
    }

    virtual double intersect(Ray &ray) = 0;
    virtual void draw() = 0;
    virtual Vector3D getNormal(const Ray &incidentRay, const Vector3D &point) = 0;

    // add proper destructor
    virtual ~Object() {}
};

class Sphere : public Object
{
public:
    Sphere(const Vector3D &center, double radius)
    {
        this->referencePoint = center;
        this->height = radius;
        this->width = radius;
        this->length = radius;
    }

    void draw() override
    {
        glColor3f(color[0], color[1], color[2]);
        glPushMatrix();
        glTranslatef(this->referencePoint.x, this->referencePoint.y, this->referencePoint.z);
        glutSolidSphere(this->length, 100, 100);
        glPopMatrix();
    }

    void setColor(double *color)
    {
        copy(color, color + 3, this->color);
        pbrMaterial.albedo = Color(color[0], color[1], color[2]);
    }

    void setShine(int shine)
    {
        this->shine = shine;
        // Map legacy shine to roughness (inverse relationship)
        pbrMaterial.roughness = 1.0 - min(1.0, shine / 100.0);
    }

    void setCoEfficients(double *coEfficients)
    {
        copy(coEfficients, coEfficients + 4, this->coEfficients);
        
        // Map legacy coefficients to PBR parameters
        pbrMaterial.ao = coEfficients[0] * 2.0;  // Ambient to AO
        pbrMaterial.roughness = 1.0 - sqrt(coEfficients[2]);  // Inverse of specular
        pbrMaterial.metallic = coEfficients[3];  // Reflection to metallic
    }

    Vector3D getNormal(const Ray &incidentRay, const Vector3D &point)
    {
        return ((this->referencePoint * -1.0) + point).normalize();
    }

    double intersect(Ray &ray) override
    {
        Vector3D q = ray.start - this->referencePoint;
        double b = 2.0 * (q.dot(ray.dir));
        double c = q.dot(q) - this->length * this->length;
        double D = b * b - 4 * c; // as a = 1 due to ray.dir.length() = 1

        double t;

        if (D < 0.0)
        {
            t = -1;
        }

        else if (D > 0.0)
        {
            double t1, t2;
            D = sqrt(D);
            t1 = (-b - D) / 2;
            t2 = (-b + D) / 2;

            if (t1 > 0)
            {
                t = t1;
            }

            else if (t2 > 0)
            {
                t = t2;
            }

            else
            {
                t = -1;
            }
        }

        else
        {
            t = -b / 2;

            if (t < 0)
            {
                t = -1;
            }
        }

        return t;
    }
};

class Triangle : public Object
{
private:
    Vector3D a, b, c;

public:
    Triangle(const Vector3D &a, const Vector3D &b, const Vector3D &c)
    {
        this->a = a;
        this->b = b;
        this->c = c;
    }

    void draw()
    {
        glColor3f(color[0], color[1], color[2]);
        glBegin(GL_TRIANGLES);
        glVertex3f(a.x, a.y, a.z);
        glVertex3f(b.x, b.y, b.z);
        glVertex3f(c.x, c.y, c.z);
        glEnd();
    }

    void setColor(double *color)
    {
        copy(color, color + 3, this->color);
    }

    void setShine(int shine)
    {
        this->shine = shine;
    }

    void setCoEfficients(double *coEfficients)
    {
        copy(coEfficients, coEfficients + 4, this->coEfficients);
    }

    Vector3D getNormal(const Ray &incidentRay, const Vector3D &point)
    {
        Vector3D normal = ((b - a).cross(c - a)).normalize();
        if (normal.dot(incidentRay.dir) > 0.0)
        {
            normal = normal * -1.0;
        }
        return normal;
    }

    double getDeterminant(double mat[3][3])
    {
        double det = 0;
        det = mat[0][0] * (mat[1][1] * mat[2][2] - mat[2][1] * mat[1][2]) - mat[0][1] * (mat[1][0] * mat[2][2] - mat[2][0] * mat[1][2]) + mat[0][2] * (mat[1][0] * mat[2][1] - mat[2][0] * mat[1][1]);
        return det;
    }

    double intersect(Ray &ray) override
    {
        double matA[3][3] = {
            {a.x - b.x, a.x - c.x, ray.dir.x},
            {a.y - b.y, a.y - c.y, ray.dir.y},
            {a.z - b.z, a.z - c.z, ray.dir.z}};

        double detA = this->getDeterminant(matA);

        assert(detA != 0.0);

        double matBeta[3][3] = {
            {a.x - ray.start.x, a.x - c.x, ray.dir.x},
            {a.y - ray.start.y, a.y - c.y, ray.dir.y},
            {a.z - ray.start.z, a.z - c.z, ray.dir.z}};

        double beta = this->getDeterminant(matBeta) / detA;

        double matGamma[3][3] = {
            {a.x - b.x, a.x - ray.start.x, ray.dir.x},
            {a.y - b.y, a.y - ray.start.y, ray.dir.y},
            {a.z - b.z, a.z - ray.start.z, ray.dir.z}};

        double gamma = this->getDeterminant(matGamma) / detA;

        double matT[3][3] = {
            {a.x - b.x, a.x - c.x, a.x - ray.start.x},
            {a.y - b.y, a.y - c.y, a.y - ray.start.y},
            {a.z - b.z, a.z - c.z, a.z - ray.start.z}};

        double t = this->getDeterminant(matT) / detA;

        if (beta + gamma < 1 && beta > 0 && gamma > 0 && t > 0)
        {
            return t;
        }
        else
        {
            return -1;
        }
    }
};

class General : public Object
{
private:
    double A, B, C, D, E, F, G, H, I, J;

public:
    // Equation: F(x,y,z) = Ax2+By2+Cz2+Dxy+Exz+Fyz+Gx+Hy+Iz+J = 0
    General(double A, double B, double C, double D, double E, double F, double G, double H, double I, double J,
            const Vector3D &referencePoint, double length, double width, double height)
    {
        this->A = A;
        this->B = B;
        this->C = C;
        this->D = D;
        this->E = E;
        this->F = F;
        this->G = G;
        this->H = H;
        this->I = I;
        this->J = J;

        this->referencePoint = referencePoint;

        this->length = length;
        this->width = width;
        this->height = height;
    }

    void draw()
    {
    }

    Vector3D getNormal(const Ray &incidentRay, const Vector3D &point)
    {
        double x = point.x;
        double y = point.y;
        double z = point.z;

        double nx = 2 * A * x + D * y + E * z + G;
        double ny = 2 * B * y + D * x + F * z + H;
        double nz = 2 * C * z + E * x + F * y + I;

        return Vector3D(nx, ny, nz).normalize();
    }

    bool isInsideBox(const Vector3D &point)
    {
        if (this->length > EPSILON)
        {
            if (point.x < this->referencePoint.x || point.x > this->referencePoint.x + this->length)
            {
                return false;
            }
        }
        if (this->width > EPSILON)
        {
            if (point.y < this->referencePoint.y || point.y > this->referencePoint.y + this->width)
            {
                return false;
            }
        }
        if (this->height > EPSILON)
        {
            if (point.z < this->referencePoint.z || point.z > this->referencePoint.z + this->height)
            {
                return false;
            }
        }
        return true;
    }

    double intersect(Ray &ray) override
    {
        double x0 = ray.start.x;
        double y0 = ray.start.y;
        double z0 = ray.start.z;

        double xd = ray.dir.x;
        double yd = ray.dir.y;
        double zd = ray.dir.z;

        double a = A * xd * xd + B * yd * yd + C * zd * zd + D * xd * yd + E * xd * zd + F * yd * zd;
        double b = 2 * (A * x0 * xd + B * y0 * yd + C * z0 * zd) + D * (x0 * yd + y0 * xd) + E * (x0 * zd + z0 * xd) + F * (y0 * zd + z0 * yd) + G * xd + H * yd + I * zd;
        double c = A * x0 * x0 + B * y0 * y0 + C * z0 * z0 + D * x0 * y0 + E * x0 * z0 + F * y0 * z0 + G * x0 + H * y0 + I * z0 + J;

        double D = b * b - 4 * a * c;

        Vector3D intersectionPoint;

        if (fabs(a) < EPSILON)
        {
            double t = -c / b;
            if (t > 0)
            {
                intersectionPoint = ray.start + ray.dir * t;
                if (this->isInsideBox(intersectionPoint))
                {
                    return t;
                }
            }
        }

        if (D < 0)
        {
            return -1;
        }

        else if (D > 0.0)
        {
            D = sqrt(D);
            double t1 = (-b - D) / (2 * a);
            double t2 = (-b + D) / (2 * a);

            if (t1 > t2)
            {
                swap(t1, t2);
            }

            if (t1 > 0)
            {
                intersectionPoint = ray.start + ray.dir * t1;
                if (this->isInsideBox(intersectionPoint))
                {
                    return t1;
                }
            }

            if (t2 > 0)
            {
                intersectionPoint = ray.start + ray.dir * t2;
                if (this->isInsideBox(intersectionPoint))
                {
                    return t2;
                }
            }
            return -1;
        }

        else
        {
            double t = -b / (2 * a);
            if (t > 0)
            {
                intersectionPoint = ray.start + ray.dir * t;
                if (this->isInsideBox(intersectionPoint))
                {
                    return t;
                }
            }
            return -1;
        }
    }
};

class Floor : public Object
{
public:
    bool useTexture;
    double textureScale;  // Controls texture repetition: <1 = more repeats, >1 = fewer repeats
    
    Floor(double floorWidth, double tileWidth)
    {
        this->referencePoint = Vector3D(-floorWidth / 2, -floorWidth / 2, 0);
        this->length = tileWidth;
        this->width = tileWidth;
        this->useTexture = false;
        this->textureScale = 1.0;  // Default: one texture per tile
    }

    void enableTexture(bool enable) {
        this->useTexture = enable;
    }
    
    void setTextureScale(double scale) {
        this->textureScale = scale;
    }

    Color sampleTexture(double u, double v) {
        if (!textureData || textureWidth <= 0 || textureHeight <= 0) {
            return Color(0.5, 0.5, 0.5); // Gray fallback
        }
        
        // Clamp u and v to [0,1]
        u = std::max(0.0, std::min(1.0, u));
        v = std::max(0.0, std::min(1.0, v));
        
        // Normalized -> pixel coords
        int pixel_x = (int)(u * (textureWidth - 1));
        int pixel_y = (int)((1.0 - v) * (textureHeight - 1)); // Flip Y
        
        // Safety clamp
        pixel_x = std::max(0, std::min(textureWidth - 1, pixel_x));
        pixel_y = std::max(0, std::min(textureHeight - 1, pixel_y));
        
        // Compute array index
        int index = (pixel_y * textureWidth + pixel_x) * textureChannels;
        int max_index = textureWidth * textureHeight * textureChannels;
        if (index < 0 || index + 2 >= max_index) {
            return Color(1.0, 0.0, 1.0); // Magenta = error
        }
        
        Color color;
        color.r = textureData[index] / 255.0;
        
        if (textureChannels >= 2) {
            color.g = textureData[index + 1] / 255.0;
        } else {
            color.g = color.r; // Grayscale
        }
        
        if (textureChannels >= 3) {
            color.b = textureData[index + 2] / 255.0;
        } else {
            color.b = color.r; // Grayscale
        }
        
        return color;
    }

    void draw()
    {
        glPushMatrix();
        glTranslatef(this->referencePoint.x, this->referencePoint.y, this->referencePoint.z);
        int tileCount = round(-this->referencePoint.x * 2 / this->length);
        for (int i = 0; i < tileCount; i++)
        {
            for (int j = 0; j < tileCount; j++)
            {
                if ((i + j) % 2 == 0)
                {
                    glColor3f(1, 1, 1);
                }
                else
                {
                    glColor3f(0, 0, 0);
                }
                glBegin(GL_QUADS);
                glVertex3f(i * this->length, j * this->length, 0);
                glVertex3f(i * this->length + this->length, j * this->length, 0);
                glVertex3f(i * this->length + this->length, j * this->length + this->length, 0);
                glVertex3f(i * this->length, j * this->length + this->length, 0);
                glEnd();
            }
        }
        glPopMatrix();
    }

    void setColor(double *color)
    {
        copy(color, color + 3, this->color);
    }

    void setShine(int shine)
    {
        this->shine = shine;
    }

    void setCoEfficients(double *coEfficients)
    {
        copy(coEfficients, coEfficients + 4, this->coEfficients);
    }

    Vector3D getNormal(const Ray &incidentRay, const Vector3D &point)
    {
        if (incidentRay.dir.z < 0.0)
        {
            return Vector3D(0, 0, 1);
        }
        return Vector3D(0, 0, -1);
    }

    double *getColor(Vector3D &point) override
    {
        if (useTexture) {
            // Calculate raw UV coordinates based on the intersection point
            double u_raw = (point.x - this->referencePoint.x);
            double v_raw = (point.y - this->referencePoint.y);
            
            // Apply texture scaling - higher scale means texture repeats less frequently
            u_raw /= (this->length * textureScale);
            v_raw /= (this->width * textureScale);
            
            // Normalize u and v for texture sampling (create repeating pattern)
            double u_norm = fmod(u_raw, 1.0);
            double v_norm = fmod(v_raw, 1.0);
            
            // Ensure positive values for fmod results
            if (u_norm < 0) u_norm += 1.0;
            if (v_norm < 0) v_norm += 1.0;
            
            return sampleTexture(u_norm, v_norm).toArray();
        } else {
            // Original checkerboard logic
            int tileX = (point.x - this->referencePoint.x) / this->length;
            int tileY = (point.y - this->referencePoint.y) / this->length;
            if ((tileX + tileY) % 2 == 0)
            {
                return whiteColor;
            }
            else
            {
                return blackColor;
            }
        }
    }

    double intersect(Ray &ray) override
    {
        assert(ray.dir.length() != 0);

        if (fabs(ray.dir.z) < EPSILON)
        {
            // The ray is parallel to the floor, no intersection
            return -1;
        }

        Vector3D normal = this->getNormal(ray, this->referencePoint);
        double t = -(normal.dot(ray.start)) / (normal.dot(ray.dir));
        if (t < 0.0)
        {
            t = -1;
        }
        else
        {
            // Compute intersection point
            Vector3D intersectionPoint = ray.start + ray.dir * t;

            // Check if the intersection point is within the floor boundaries
            int tileCount = round(-this->referencePoint.x * 2 / this->length);
            double minX = this->referencePoint.x;
            double minY = this->referencePoint.y;
            double maxX = minX + tileCount * this->length;
            double maxY = minY + tileCount * this->length;

            if (intersectionPoint.x >= minX && intersectionPoint.x <= maxX &&
                intersectionPoint.y >= minY && intersectionPoint.y <= maxY)
            {
                return t;
            }
            else
            {
                t = -1;
            }
        }
        return t;
    }
};