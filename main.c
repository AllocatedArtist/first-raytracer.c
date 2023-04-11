#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#include <stdlib.h>


/* ====================================================================================================================
    pretty much a c-ified version of c++ scratchapixel basic raytracer, with my own notes for learning

    why?

    this is a way for me to gain more understanding in the raytracing algorithm discussed. and to see something cool :)

    inspiration: https://github.com/scratchapixel/code/blob/main/introduction-to-ray-tracing/raytracer.cpp
    ====================================================================================================================
*/ 

#include <assert.h>

#define PI 3.141592653589793
#define PRINT(x) printf("%.4f\n", x)

typedef struct {
    float x, y, z;
} vector3f_t;

vector3f_t vector3f_empty() { return (vector3f_t) { 0.0, 0.0, 0.0 }; }
vector3f_t vector3f_new_one_val(float val) { return (vector3f_t) { val, val, val }; }
vector3f_t vector3f_new(float x, float y, float z) { return (vector3f_t) { x, y, z }; }
float vector3f_length2(vector3f_t vec) { return (vec.x*vec.x) + (vec.y*vec.y) + (vec.z*vec.z); } //length squared ^2
float vector3f_length_sqr(vector3f_t vec) { return sqrtf(vector3f_length2(vec)); } //square root of length
vector3f_t vector3f_normalize(vector3f_t vec) { float m = vector3f_length_sqr(vec); return m > 0 ? (vector3f_t) { vec.x / m, vec.y / m, vec.z / m } : vec; }
vector3f_t vector3f_multiply_f(vector3f_t vec, float multiplier) {  return (vector3f_t) { vec.x * multiplier, vec.y * multiplier, vec.z * multiplier }; }
vector3f_t vector3f_multiply_vec(vector3f_t a, vector3f_t b) { return (vector3f_t) { a.x * b.x, a.y * b.y, a.z * b.z }; }
float vector3f_dot_product(vector3f_t a, vector3f_t b) { return (a.x * b.x) + (a.y * b.y) + (a.z * b.z); }
vector3f_t vector3f_add(vector3f_t a, vector3f_t b) { return (vector3f_t) { a.x + b.x, a.y + b.y, a.z + b.z }; }
vector3f_t vector3f_subtract(vector3f_t a, vector3f_t b) { return (vector3f_t) { a.x - b.x, a.y - b.y, a.z - b.z }; }

void vector3f_print(vector3f_t vec) { printf("(%.4f, %.4f, %.4f)\n", vec.x, vec.y, vec.z ); }
bool float_is_equal(float x, float y) { if (fabs(x - y) < 0.01) return true; else return false; } 
void vector3f_unit_tests(void) { 
    vector3f_t v1 = vector3f_empty();
    assert(v1.x == 0.0 && v1.y == 0.0 && v1.z == 0.0);

    vector3f_t v2 = vector3f_new_one_val(69.0);
    assert(v2.x == 69.0 && v2.y == 69.0 && v2.z == 69.0);

    vector3f_t v3 = vector3f_normalize(vector3f_new(21.0, 420.0, 69.0));
    assert(float_is_equal(v3.x, (float)0.049279) && float_is_equal(v3.y, (float)0.98557) && float_is_equal(v3.z, (float)0.161916));

    vector3f_t v4 = vector3f_multiply_f(vector3f_new_one_val(8.0), 7.0);
    assert(v4.x == 56.0 && v4.y == 56.0 && v4.z == 56.0);

    vector3f_t v5 = vector3f_multiply_vec(vector3f_new_one_val(9.0), vector3f_new_one_val(3.0));
    assert(v5.x == 27.0 && v5.y == 27.0 && v5.z == 27.0);

    float dot_product = vector3f_dot_product(vector3f_new(1.0, -46.0, 9.0), vector3f_new(-8.0, 22.0, 13.0));
    assert(dot_product == -903.0);

    vector3f_t v6 = vector3f_add(vector3f_new(39.0, 12.0, 2.0), vector3f_new_one_val(55.0));
    assert(v6.x == 94.0 && v6.y == 67.0 && v6.z == 57.0);

    vector3f_t v7 = vector3f_subtract(vector3f_new(39.0, 12.0, 2.0), vector3f_new_one_val(55.0));
    assert(v7.x == -16.0 && v7.y == -43.0 && v7.z == -53.0);
}

typedef struct {
    vector3f_t center; //sphere position
    float radius1, radius2; //radius & radius squared ^2
    vector3f_t surface_color, emission_color;
    float transparency, reflection; 
} sphere_t;

sphere_t sphere_new(vector3f_t center, float radius, vector3f_t surface_c, vector3f_t emission_c, float transparency, float reflection) {
    return (sphere_t) { center, radius, powf(radius, 2), surface_c, emission_c, transparency, reflection };
}

//ray_dir should be normalized when passing into the function **using geometric solution**
/*

    a ray can be defined as origin + (t * dir)
    if t is equal to 0, the end of the ray will be the same as the origin
    if t is greater than 0, the end of the ray will be in front of the origin
    if t is less than 0, the end of the ray will be behind the origin

    t0 and t1 (float values) are the distance between the ray origin and the points of intersection (if any) on the sphere
    t0 is the closest intersection and t1 is the farthest

    take P and P' where P is the closest point of intersection (this time a vector3) and P' is the farthest point of intersection
    once we have t0 and t1, we can find the exact location of the intersections with:
    
    P  (closest point as vector3) = origin + (t0 * dir)
    P' (farthest point as vector3) = origin + (t1 * dir)

    so t0 and t1 are important.

    now in order to get t0 and t1 we need to find L, tca, and thc. 

    L is the vector between the ray origin and the sphere center. basically it is a vector that points from 
    ray origin to the sphere center.

    tca is the dot product between L and D (ray direction). this is the same as projecting
    D onto L (scalar projection). what this means is that we get the length
    of a segment between the ray direction and the center of the sphere. 
    
    if tca is negative then that means the ray direction is not facing towards the sphere at all!
    so that means no intersection is possible. to be more specific, the ray faces away from the sphere.
    the sphere is behind the ray or vice versa :)

    simply: x > 0: facing towards the sphere. x == 0 facing exactly towards the sphere. x < 0: facing away from the sphere

    to find thc we will need to find d (not D).

    basically, d, tca, and L form a right triangle.
    we have tca and L, so to find d we use the pythagorean theoreom. 
    L is the hypotenuse, so the equation looks like this d^2 + tca^2 = L^2
    Using basic algebra we can: d^2 = L^2 - tca^2
    d = sqr(L^2 - tca^2)

    wait! but L is a vector, and we are trying to subtract that by a float? well we just 
    find the magnitude/length squared of the vector. so more accurately it is:

    d = sqr(length squared of L - tca^2)
    
    if d is greater than the sphere's radius, then that means the ray isn't intersecting.
    the ray faces the same direction as the sphere, but it shoots outside and therefore
    doesn't intersect. if radius is the size of the sphere, then we would have to be 
    inside, or within the sphere's radius (less than radius) in order to be intersecting 

    to optimize, we could scrap the square root and just take (length of L squared - tca ^ 2).
    then we just compare to radius squared. this avoids expensive square root calculation when
    we don't need it.

    now to get thc we recognize that d and radius form a right triangle with thc.
    so the equation is again thc + d = r
    thc = r - d

    now for this one we would like to take the square root.
    so it is the square root of radius squared - d squared

    we still end up using a sqrt anyways... :(

    anyways t0 = tca - thc
    t1 = tca + thc

    we can think of tca as the length from the origin to how much it "penetrates into the sphere" :D
    and thc is the rest of the sphere that comes after tca

    so obviously we can get P' by adding tca + thc (the full length)
    and then P by subtracting tca - thc (first contact)

*/
bool sphere_v_ray(sphere_t sphere, vector3f_t origin, vector3f_t ray_dir, float* t0, float* t1) {
    vector3f_t l = vector3f_subtract(sphere.center, origin);

    float tca = vector3f_dot_product(l, ray_dir);

    if (tca < 0.0) { return false; }

    float d2 = vector3f_length2(l) - (tca * tca); //d2 because it is squared and we never took square root

    if (d2 > sphere.radius2) { return false; }

    float thc = sqrtf(sphere.radius2 - d2);

    *t0 = tca - thc;
    *t1 = tca + thc;

    return true; 
}

//Amount of reflections ooooooh (let's burn computer D:)
#define MAX_RAY_DEPTH 5

//tweak for different effects
#define FRESNEL_MIX_VALUE 0.1

//slight offset in reflection/refraction/transmission -- also avoid hitting same sphere again when tracing ...i think
#define BIAS 1e-4

//mix function in scratchapixel demo. but it's basically a lerp... soo 
float lerp(float a, float b, float t) { return b * t + a * (1.0 - t); }
float max(float a, float b) { return a > b ? a : b; }
float min(float a, float b) { return a < b ? a : b; }

vector3f_t trace(sphere_t* spheres, size_t length, vector3f_t origin, vector3f_t ray_dir, int depth) {

    //assert(vector3f_length_sqr(ray_dir) == 1.0000); //force this to be normalized >:)

    //tnear is the closest object. at the moment we haven't hit anything so tnear is infinite
    float tnear = INFINITY;
    const sphere_t* sphere = NULL; //closest sphere

    for (int i = 0; i < length; ++i) {
        float t0 = INFINITY, t1 = INFINITY;
        if (sphere_v_ray(spheres[i], origin, ray_dir, &t0, &t1)) {
            if (t0 < 0.0) t0 = t1; //the ray origin is inside the sphere so closest ray intersection is now t1
            if (t0 < tnear) { //found new closest distance and object
                tnear = t0;
                sphere = spheres + i;
            }
        }
    }

    if (sphere == NULL) return vector3f_new_one_val(2.0); //no intersection return black or background color

    vector3f_t surface_color = vector3f_empty();
    vector3f_t p_hit = vector3f_add(origin, vector3f_multiply_f(ray_dir, tnear)); //remember origin + (t0 * dir) = P (closest/first point of intersection)?
    vector3f_t p_normal = vector3f_normalize(vector3f_subtract(p_hit, sphere->center)); //get direction from center to P and then normalize. 

    bool inside = false;
    if (vector3f_dot_product(ray_dir, p_normal) > 0.0) {
        p_normal = vector3f_multiply_f(p_normal, -1.0); //reverse normal if the normal and view face similar directions
        inside = true;
    }

    //cool effects happen
    if ((sphere->transparency > 0.0 || sphere->reflection > 0.0) && depth < MAX_RAY_DEPTH) {

        //cosine of incoming ray and normal? probably the angle at which the ray hits the sphere. 
        //which distorts slightly as the ray enters the refractive/reflective sphere and as it
        //exits. of course this is why we want the ray direction to initially face opposite to the normal 
        //this is why we reverse the normal direction when the ray direction happens to be the same because the ray is in the sphere
        float facing_ratio = -vector3f_dot_product(ray_dir, p_normal);


        //fresnel equation gives us mix of refraction to reflection
        float fresnel_effect = lerp(pow(1.0 - facing_ratio, 3.0), 1.0, FRESNEL_MIX_VALUE);

        //reflection formula for a surface with a normal is ray_dir - normal * dot(ray_dir, normal) * 2
        //no need to normalize because we assume all values have already been normalized. which they have hopefully :)
        vector3f_t reflect_dir = vector3f_subtract(ray_dir, vector3f_multiply_f(p_normal, vector3f_dot_product(ray_dir, p_normal) * 2.0)); 

        //recursion because there can be multiple reflections of multiple reflective objects. also we set the origin slightly above the 
        //initial point of contact because of precision errors and stuff. also because it would be bad if we hit the sphere we already traced
        vector3f_t reflection = trace(spheres, length, vector3f_add(p_hit, vector3f_multiply_f(p_normal, BIAS)), reflect_dir, depth + 1);

        vector3f_t refraction = vector3f_empty();
        if (sphere->transparency == 1.0) {
            //ior = index of refraction. eta is a ratio. if we are inside the sphere then eta is ior which is the 
            //index of refraction when inside. outside has a different index of refraction so we take the ratio
            //of the outside refraction to the inside refraction. according to scratchapixel, the refraction of air
            //is very close to 1 whereas glass is very close to 1.5. we can tweak index of refraction for different
            //effects.
            float ior = 1.1, eta = (inside) ? ior : 1.0 / ior;
            float cosi = -vector3f_dot_product(p_normal, ray_dir); //cos inverse

            //reflection coefficient. determines how much light reflected from a surface
            //contributes to final color. smaller eta value means more reflection whereas
            //eta over 1 value means no reflection at all. TODO: UNDERSTAND HOW/WHAT/WHY THE HECK K AND COSI MEAN/CALCULATE/WHATEVER
            float k = 1.0 - eta * eta * (1.0 - cosi * cosi);

            vector3f_t a = vector3f_multiply_f(ray_dir, eta);
            vector3f_t b = vector3f_multiply_f(p_normal, eta * cosi - sqrtf(k));

            vector3f_t refraction_dir = vector3f_add(a, b);
            refraction_dir = vector3f_normalize(refraction_dir);

            //new direction after passing through intersection of object
            refraction = trace(spheres, length, vector3f_subtract(p_hit, vector3f_multiply_f(p_normal, BIAS)), refraction_dir, depth + 1);
        }

        vector3f_t x = vector3f_multiply_f(reflection, fresnel_effect);
        vector3f_t y = vector3f_multiply_f(refraction, sphere->transparency * (1.0 - fresnel_effect));
        vector3f_t z = vector3f_add(x, y);
        surface_color = vector3f_multiply_vec(z, sphere->surface_color);
    } else {
        //diffuse object with no reflection/refractive properties 
        for (int i = 0; i < length; ++i) {
            if (spheres[i].emission_color.x > 0.0) {
                vector3f_t transmission = vector3f_new(1.0, 1.0, 1.0); //light color
                vector3f_t light_dir = vector3f_normalize(vector3f_subtract(spheres[i].center, p_hit));

                for (int j = 0; j < length; ++j) {
                    if (i != j) {
                        float t0, t1;
                        //shadows because sphere[i] is being blocked by sphere[j] from the light
                        if (sphere_v_ray(spheres[j], vector3f_add(p_hit, vector3f_multiply_f(p_normal, BIAS)), light_dir, &t0, &t1)) {
                            transmission = vector3f_new(0.0, 0.0, 0.0);
                            break;
                        }
                    }
                }
                
                //oh boy...
                vector3f_t x = vector3f_multiply_vec(sphere->surface_color, transmission);
                float y = max(0.0, vector3f_dot_product(p_normal, light_dir));
                vector3f_t z = vector3f_multiply_f(x, y);

                surface_color = vector3f_add(surface_color, vector3f_multiply_vec(z, spheres[i].emission_color));
                
            }
        }
    }

    return vector3f_add(surface_color, sphere->emission_color);
}

void render(sphere_t* spheres, size_t length) {
    unsigned int width = 800, height = 640;
    vector3f_t* image = malloc((width* height) * sizeof(vector3f_t));

    if (image == NULL) {
        printf("COULD NOT ALLOCATE IMAGE\n");
        exit(-1);
    }

    vector3f_t* pixel = image;

    float inv_width = 1.0 / (float)width, inv_height = 1.0 / (float)height;
    float fov = 40.0, aspect_ratio = width / (float)height;
    float angle = tan(PI * 0.5 * fov / 180.0);

    //TODO UNDERSTAND THIS:
    for (unsigned int y = 0; y < height; ++y) {
        for (unsigned int x = 0; x < width; ++x, ++pixel) {
            float xx = (2.0 * ((x + 0.5) * inv_width) - 1.0) * angle * aspect_ratio;
            float yy = (1.0 - 2.0 * ((y + 0.5) * inv_height)) * angle;
            vector3f_t ray_dir = vector3f_normalize(vector3f_new(xx, yy, -1.0));
            *pixel = trace(spheres, length, vector3f_empty(), ray_dir, 0);
        }
    }

    FILE* fp = fopen("image.ppm", "wb");
    if (fp == NULL) {
        printf("ERROR: COULD NOT OPEN FILE.\n");
        exit(-1);
    }
    
    fprintf(fp, "P6\n%i %i\n255\n", width, height);

    for (unsigned int i = 0; i < width * height; ++i) {

        unsigned char r = (unsigned char)(min(1.0, image[i].x) * 255.0);
        unsigned char g = (unsigned char)(min(1.0, image[i].y) * 255.0);
        unsigned char b = (unsigned char)(min(1.0, image[i].z) * 255.0);

        //fwrite(&r, sizeof(unsigned char), 1, fp);
         //fwrite(&g, sizeof(unsigned char), 1, fp); //alternative
        //fwrite(&b, sizeof(unsigned char), 1, fp);
         fprintf(fp, "%c%c%c", r, g, b); //DON'T SPACE %c out OR ELSE IT DOESN'T WRITE PROPERLY!!!
    }


    fclose(fp);

    free(image);
}

int main(void) {

    srand(3);
    
    sphere_t spheres[6];

    spheres[0] = sphere_new(vector3f_new(-4.0, 0.0, -20.0), 3.0, vector3f_new(0.70, 0.32, 0.36), vector3f_empty(), 1.0, 0.0);
    spheres[1] = sphere_new(vector3f_new(5.0, 3.0, -60.0), 6.0, vector3f_new(0.50, 0.32, 1.00), vector3f_empty(), 1.0, 0.0);
    spheres[2] = sphere_new(vector3f_new(0, -1030.0, -100.0), 1000.0, vector3f_new(0.3, 0.3, 0.3), vector3f_empty(), 0.0, 1.0);
    spheres[3] = sphere_new(vector3f_new(3.0, 0.0, -10.0), 0.8, vector3f_new(1.0, 1.0, 0.0), vector3f_empty(), 0.0, 1.0);
    spheres[4] = sphere_new(vector3f_new(0.0, 0.0, -70.0), 10.0, vector3f_new(0.43, 1.0, 0.32), vector3f_new(0.0, 0.0, 0.0), 1.0, 1.0);
    spheres[5] = sphere_new(vector3f_new(0.0, 20.0, -30.0), 3.0, vector3f_new_one_val(0.0), vector3f_new(0.5, 0.5, 0.5), 1.0, 0.0);

    size_t length = sizeof(spheres) / sizeof(spheres[0]);

    //sphere_t example_spheres[6];
   //example_spheres[0] = sphere_new(vector3f_new(0.0, -10004.0, -20.0), 10000.0, vector3f_new(0.20, 0.20, 0.20), vector3f_empty(), 0.0, 0.0);
   //example_spheres[1] = sphere_new(vector3f_new(0.0, 0.0, -20.0), 4.0, vector3f_new(1.00, 0.32, 0.36), vector3f_empty(), 1.0, 0.0);
   //example_spheres[2] = sphere_new(vector3f_new(5.0, -1.0, -15.0), 2.0, vector3f_new(0.9, 0.76, 0.46), vector3f_empty(), 1.0, 0.0);
   //example_spheres[3] = sphere_new(vector3f_new(5.0, 0.0, -25.0), 3.0, vector3f_new(0.65, 0.77, 0.97), vector3f_empty(), 1.0, 0.0);
   //example_spheres[4] = sphere_new(vector3f_new(-5.5, 0.0, -15.0), 3.0, vector3f_new(0.9, 0.9, 0.9), vector3f_empty(), 1.0, 0.0);
   //example_spheres[5] = sphere_new(vector3f_new(0.0, 20.0, -30.0), 3.0, vector3f_new_one_val(0.0), vector3f_new(3.0, 3.0, 3.0), 0.0, 0.0);

    render(spheres, length);
    

    return 0;
}