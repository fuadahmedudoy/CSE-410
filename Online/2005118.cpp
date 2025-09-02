/**
 * OpenGL 3D Drawing Demo
 *
 * This program demonstrates basic 3D rendering with OpenGL and GLUT including:
 * - Camera positioning with gluLookAt
 * - Drawing 3D shapes (cube and pyramid)
 * - Keyboard navigation for camera control
 * - Perspective projection
 * - Object toggling
 */

// --- Includes ---
// Standard Headers
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// OpenGL / GLUT Headers
#ifdef __APPLE__
#include <GLUT/glut.h> // Use GLUT framework on macOS
#else
#include <GL/glut.h> // Use standard GLUT location on Linux/Windows
#endif

#define PI 3.14159265358979323846 // Define PI constant

// --- Global Variables ---
// Camera position and orientation
GLfloat eyex = 4, eyey = 4, eyez = 4;          // Camera position coordinates
GLfloat centerx = 0, centery = 0, centerz = 0; // Look-at point coordinates
GLfloat upx = 0, upy = 1, upz = 0;             // Up vector coordinates
GLfloat spherex,spherey,spherez; // Sphere coordinates
GLfloat velocityx=2,velocityy=4,velocityz=2; // Sphere velocity
GLfloat rotationx=0,rotationy=0,rotationz=0; // Sphere rotation angles
GLfloat velocity=5.0f,yawAngle=0.0f; // Sphere velocity and yaw angle
GLfloat rotationAxisx=0,rotationAxisy=0,rotationAxisz=0,rotationAngle=0,bladeang=0,bladespeed=1,fanang=0; // Sphere rotation axis
// Object visibility flags
bool isAxes = false;     // Toggle for coordinate axes
bool isCube = true;    // Toggle for cube
bool isPyramid = false; // Toggle for pyramid
bool isMovingSphere = false,isReset=true; // Toggle for moving sphere
bool showVector = false,isYaw=false; // Toggle for vector display

// --- Function Declarations ---
void initGL();
void display();
void reshapeListener(GLsizei width, GLsizei height);
void keyboardListener(unsigned char key, int x, int y);
void specialKeyListener(int key, int x, int y);
void drawAxes();
void drawCube();
void drawPyramid();
void spherePosition();
void drawVelocityVector();
void calculateVelocity();
/**
 * Initialize OpenGL settings
 * Sets up background color and enables depth testing
 */
void initGL()
{
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f); // Black background
    glEnable(GL_DEPTH_TEST);              // Enable depth testing for z-culling
}
void drawAxes()
{
    glLineWidth(3); // Set line thickness

    glBegin(GL_LINES);

    // X axis (red)
    glColor3f(1, 1, 1);
    glVertex3f(-100, 0, 0);
    glVertex3f(100, 0, 0);

    // Y axis (green)
    glColor3f(1, 1, 1);
    glVertex3f(0, -100, 0);
    glVertex3f(0, 100, 0);

    // Z axis (blue)
    glColor3f(1, 1, 1);
    glVertex3f(0, 0, -100);
    glVertex3f(0, 0, 100);

    glEnd();
}
void drawCylinder(float baseRadius, float topRadius, float height)
{
    GLUquadric *quad = gluNewQuadric();
    glColor3f(0.675f,0.455f,0.910f); // Gray hub
    gluCylinder(quad, baseRadius, topRadius, height, 20, 20);
    gluDeleteQuadric(quad);
}
void drawBlade(){
    glPushMatrix();
    glTranslatef(0,0,0.75);
    glLineWidth(8); // Set line thickness
    glBegin(GL_LINES);
    glColor3f(0.53f, 0.81f, 0.92f);
    glVertex3f(0,0,0);
    glVertex3f(0,2,0);
    glEnd();
    glPopMatrix();

}
void drawRotatingBlade(){
   for(int i=0;i<6;i++){
    glPushMatrix();
    glRotatef(i*60,0,0,1);
    drawBlade();
    glPopMatrix();
   }
}
void drawBox()
{
    glPushMatrix();
    glTranslatef(0,1.5,0);
    glPushMatrix();
        float length=1,width=1,height=1;
        glBegin(GL_QUADS);
        glColor3f(0.675f,0.455f,0.910f);
        glVertex3f( length,  width, -height);
        glVertex3f(-length,  width, -height);
        glVertex3f(-length,  width,  height);
        glVertex3f( length,  width,  height);
        glEnd();
    glPopMatrix();
    glPushMatrix();
        glTranslatef(-1,-0.4,0);
        glRotatef(45,0,0,1);
        length=1,width=1,height=1;
        glBegin(GL_QUADS);
        glColor3f(0.941f,0.322f,0.322f);
        glVertex3f( length,  width, -height);
        glVertex3f(-length,  width, -height);
        glVertex3f(-length,  width,  height);
        glVertex3f( length,  width,  height);
        glEnd();
    glPopMatrix();
    glTranslatef(1,-0.4,0);
    glPushMatrix();
        glRotatef(-45,0,0,1);
        length=1,width=1,height=1;
        glBegin(GL_QUADS);
        glColor3f(0.941f,0.322f,0.322f);
        glVertex3f( length,  width, -height);
        glVertex3f(-length,  width, -height);
        glVertex3f(-length,  width,  height);
        glVertex3f( length,  width,  height);
        glEnd();
    glPopMatrix();
    glTranslatef(0.42,-1,0);
    glPushMatrix();
        glRotatef(-90,0,0,1);
        length=1,width=1,height=1;
        glBegin(GL_QUADS);
        glColor3f(0.675f,0.455f,0.910f);
        glVertex3f( length,  width, -height);
        glVertex3f(-length,  width, -height);
        glVertex3f(-length,  width,  height);
        glVertex3f( length,  width,  height);
        glEnd();
    glPopMatrix();
    glPushMatrix();
        glTranslatef(-4.8,0,0);
        glRotatef(-90,0,0,1);
        length=1,width=1,height=1;
        glBegin(GL_QUADS);
        glColor3f(0.675f,0.455f,0.910f);
        glVertex3f( length,  width, -height);
        glVertex3f(-length,  width, -height);
        glVertex3f(-length,  width,  height);
        glVertex3f( length,  width,  height);
        glEnd();
    glPopMatrix();
    glTranslatef(-0.4,-1,0);
    glPushMatrix();
        glRotatef(-135,0,0,1);
        length=1,width=1,height=1;
        glBegin(GL_QUADS);
        glColor3f(0.941f,0.322f,0.322f);
        glVertex3f( length,  width, -height);
        glVertex3f(-length,  width, -height);
        glVertex3f(-length,  width,  height);
        glVertex3f( length,  width,  height);
        glEnd();
    glPopMatrix();
    glPushMatrix();
        glTranslatef(-2,0,0);
        glRotatef(135,0,0,1);
        length=1,width=1,height=1;
        glBegin(GL_QUADS);
        glColor3f(0.941f,0.322f,0.322f);
        glVertex3f( length,  width, -height);
        glVertex3f(-length,  width, -height);
        glVertex3f(-length,  width,  height);
        glVertex3f( length,  width,  height);
        glEnd();
    glPopMatrix();
    glTranslatef(-0.95,-0.4,0);
    glPushMatrix();
        glRotatef(-180,0,0,1);
        length=1,width=1,height=1;
        glBegin(GL_QUADS);
        glColor3f(0.675f,0.455f,0.910f);
        glVertex3f( length,  width, -height);
        glVertex3f(-length,  width, -height);
        glVertex3f(-length,  width,  height);
        glVertex3f( length,  width,  height);
        glEnd();
    glPopMatrix();
    glPopMatrix();
}
/**
 * Main display function
 * Sets up the camera and renders visible objects
 */
void display()
{
    // Clear color and depth buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // Set up the model-view matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // glRotatef(yawAngle, upx, upy, upz);

    // Position camera using the eye, center and up vectors
    gluLookAt(eyex, eyey, eyez,          // Camera position
              centerx, centery, centerz, // Look-at point
              upx, upy, upz);            // Up vector
     // Apply yaw rotation
    // Draw objects based on visibility flags
    // if (isCube)
    //     drawCube();
   
    drawAxes();
    glPushMatrix();
        glRotatef(fanang,1,0,0);
        drawBox();
        drawCylinder(0.1f, 0.1f, 0.75f);
        glRotatef(bladeang,0,0,1);
        drawRotatingBlade();
        
    glPopMatrix();
        
    // // Swap buffers (double buffering)
    // if (showVector)
    //     drawVelocityVector();
    
    glutSwapBuffers();
}

/**
 * Window reshape callback
 * Handles window resizing and maintains aspect ratio
 */
void reshapeListener(GLsizei width, GLsizei height)
{
    // Prevent division by zero
    if (height == 0)
        height = 1;

    // Calculate aspect ratio
    GLfloat aspect = (GLfloat)width / (GLfloat)height;

    // Set viewport to cover entire window
    glViewport(0, 0, width, height);

    // Set up perspective projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    // 45-degree field of view, aspect ratio, near and far clipping planes
    gluPerspective(45.0f, aspect, 0.1f, 100.0f);
}
void yaw(float angle) {
    float dirX = centerx - eyex;
    float dirY = centery - eyey;
    float dirZ = centerz - eyez;
    float dirLen = sqrt(dirX*dirX + dirY*dirY + dirZ*dirZ);
    dirX /= dirLen; // Normalize direction vector   
    dirY /= dirLen;
    dirZ /= dirLen;
    // Use current up vector as axis of rotation
    float upLen = sqrt(upx*upx + upy*upy + upz*upz);
    float ux = upx / upLen;
    float uy = upy / upLen;
    float uz = upz / upLen;
    
    float rad = (angle * PI) / 180.0f;
    float cosA = cos(rad);
    float sinA = sin(rad);
    
    // Rodrigues' rotation formula
    float dot = dirX*ux + dirY*uy + dirZ*uz;
    
    float newDirX = dirX * cosA + (uy * dirZ - uz * dirY) * sinA ;
    float newDirY = dirY * cosA + (uz * dirX - ux * dirZ) * sinA ;
    float newDirZ = dirZ * cosA + (ux * dirY - uy * dirX) * sinA ;

    
    // Update look-at point
    centerx = eyex + newDirX;
    centery = eyey + newDirY;
    centerz = eyez + newDirZ;
    // isYaw=true;
    // yawAngle+=angle;    
}

void pitch(float angle) {
    float dirX = centerx - eyex;
    float dirY = centery - eyey;
    float dirZ = centerz - eyez;

    // Right = dir × up
    float rightX = dirY * upz - dirZ * upy;
    float rightY = dirZ * upx - dirX * upz;
    float rightZ = dirX * upy - dirY * upx;

    // Normalize right vector
    float len = sqrt(rightX * rightX + rightY * rightY + rightZ * rightZ);
    rightX /= len;
    rightY /= len;
    rightZ /= len;

    float rad = angle * PI / 180.0f;
    float cosA = cos(rad);
    float sinA = sin(rad);

    // Rotate dir using Rodrigues' formula
    float newDirX = dirX * cosA + (rightY * dirZ - rightZ * dirY) * sinA;
    float newDirY = dirY * cosA + (rightZ * dirX - rightX * dirZ) * sinA ;
    float newDirZ = dirZ * cosA + (rightX * dirY - rightY * dirX) * sinA ;

    centerx = eyex + newDirX;
    centery = eyey + newDirY;
    centerz = eyez + newDirZ;

    // Also rotate up vector around the same right axis
    float newUpX = upx * cosA + (rightY * upz - rightZ * upy) * sinA ;
    float newUpY = upy * cosA + (rightZ * upx - rightX * upz) * sinA ;
    float newUpZ = upz * cosA + (rightX * upy - rightY * upx) * sinA ;

    upx = newUpX;
    upy = newUpY;
    upz = newUpZ;
    
}


void tilt(float angle){
// 1. Direction vector
    float dirX = centerx - eyex;
    float dirY = centery - eyey;
    float dirZ = centerz - eyez;

    // Normalize dir
    float len = sqrt(dirX*dirX + dirY*dirY + dirZ*dirZ);
    dirX /= len; dirY /= len; dirZ /= len;

    // Convert angle to radians
    float rad = angle * PI / 180.0f;
    float cosA = cos(rad);
    float sinA = sin(rad);

    // Rodrigues' rotation formula
    float newupx = upx * cosA + (dirY * upz - dirZ * upy) * sinA;
    float newupy = upy * cosA + (dirZ * upx - dirX * upz) * sinA;
    float newupz = upz * cosA + (dirX * upy - dirY * upx) * sinA;

    upx = newupx;
    upy = newupy;
    upz = newupz;
}
/**
 * Keyboard input handler for standard keys
 * Manages camera position, object visibility, and program exit
 */
void keyboardListener(unsigned char key, int x, int y)
{
    float v = 0.1; // Movement increment

    switch (key)
    {
    // --- Camera Position Controls (eye coordinates) ---
    case '1':
        yaw(2);
        break; // Move eye right
    case '2':
        yaw(-2);
        break; // Move eye left
    case '3':
        pitch(2);
        break; // Move eye up
    case '4':
        pitch(-2);
        break; // Move eye down
    case '5':
        tilt(2);
        break; // Move eye forward
    case '6':
        tilt(-2);
        break; // Move eye backward
    case 'w':
        // eyex += upx*v;
        // eyey += upy*v;
        // eyez += upz*v;
        bladespeed+=2;
        break; // Move look-at point left
    case 's':
        // eyex -= upx*v;
        // eyey -= upy*v;
        // eyez -= upz*v;
        bladespeed-=2;
        break;
    case 'a':
        // eyex += upx*v;
        // eyey += upy*v;
        // eyez += upz*v;
        fanang+=5;
        break; // Move look-at point left
    case 'd':
        // eyex -= upx*v;
        // eyey -= upy*v;
        // eyez -= upz*v;
        fanang-=5;
        break;
    case 'r':
        if(!isMovingSphere){
            isReset=true;
            velocity=2.0f;
            spherePosition();
        }
        
        break;
    case 'v':
        if(!isMovingSphere){
            showVector = !showVector;
            
        }
        break;

    // --- Object Visibility Toggles ---
    case 'c':
        isCube = !isCube;
        break; // Toggle cube
    case '+':
        if(!isMovingSphere && isReset){
            velocityy += 0.5f;
           // printf("+  & velocityy: %f\n",velocityy);
            //calculateVelocity();
             // Increase velocity
        }
        break;
        
    case '-':
        if(!isMovingSphere && isReset){
           // if (velocityy > 0.5f) // Prevent negative velocity
                velocityy -= 0.5f;
               // printf("+  & velocityy: %f\n",velocityy);
                //calculateVelocity();
        }   // Decrease velocity}}   
        break; // Decrease velocity

    // --- Program Control ---
    case 32: // Space key: toggle sphere movement
        isMovingSphere = !isMovingSphere;
        isReset=false;
        //printf("Sphere movement %s\n", isMovingSphere ? "enabled" : "disabled");
        break;
    }

    glutPostRedisplay(); // Request a screen refresh
}

/**
 * Special key input handler (arrow keys, function keys)
 * Provides camera orbit functionality
 */
void moveLeftRight(int sign){
    // 1. Direction vector
    float dirX = centerx - eyex;
    float dirY = centery - eyey;
    float dirZ = centerz - eyez;

    // 2. Normalize
    float len = sqrt(dirX*dirX + dirY*dirY + dirZ*dirZ);
    dirX /= len; dirY /= len; dirZ /= len;

    // 3. Compute right vector = dir × up
    float rightX = dirY * upz - dirZ * upy;
    float rightY = dirZ * upx - dirX * upz;
    float rightZ = dirX * upy - dirY * upx;

    // 4. Normalize right
    len = sqrt(rightX*rightX + rightY*rightY + rightZ*rightZ);
    rightX /= len; rightY /= len; rightZ /= len;

    // 5. Move left = negative of right
    double v = 0.25;
    eyex += sign*(rightX * v);
    eyey += sign*(rightY * v);
    eyez += sign*(rightZ * v);

    centerx += sign*(rightX * v);
    centery += sign*(rightY * v);
    centerz += sign*(rightZ * v);

}
void forwardBackward(int sign){
    // 1. Direction vector
    float dirX = centerx - eyex;
    float dirY = centery - eyey;
    float dirZ = centerz - eyez;

    // 2. Normalize
    float len = sqrt(dirX*dirX + dirY*dirY + dirZ*dirZ);
    dirX /= len; dirY /= len; dirZ /= len;

    // 3. Move forward = direction
    double v = 0.2;
    eyex += sign*(dirX * v);
    eyey += sign*(dirY * v);
    eyez += sign*(dirZ * v);

    centerx += sign*(dirX * v);
    centery += sign*(dirY * v);
    centerz += sign*(dirZ * v);

}   
void specialKeyListener(int key, int x, int y)
{
    double v = 0.25; // Movement increment

    // Calculate view direction vector
    double lx = centerx - eyex;
    double lz = centerz - eyez;
    double s;

    switch (key)
    {
    case GLUT_KEY_LEFT:
        // Orbit camera left around the look-at point
        moveLeftRight(-1);
        break;

    case GLUT_KEY_RIGHT:
        // Orbit camera right around the look-at point
        // eyex += v * (-upy * lz);
        // eyez += v * (lx * upy);
        // s = sqrt(eyex * eyex + eyez * eyez) / (4 * sqrt(2));
        // eyex /= s;
        // eyez /= s;
        moveLeftRight(1);
        break;

    case GLUT_KEY_UP:
        // Move camera up
        forwardBackward(1);
        break;

    case GLUT_KEY_DOWN:
        // Move camera down
        forwardBackward(-1);
        break;
        
    case GLUT_KEY_PAGE_UP:
        // Move camera up
        eyey += upy*v;
        eyex += upx*v;  
        eyez += upz*v;

        centerx += upx*v;
        centery += upy*v;
        centerz += upz*v;
        break; // Move camera up  
    case GLUT_KEY_PAGE_DOWN:
        // Move camera down
        eyey -= upy*v;
        eyex -= upx*v;  
        eyez -= upz*v;

        centerx -= upx*v;
        centery -= upy*v;
        centerz -= upz*v;
        break; // Move camera down  
    }


    glutPostRedisplay(); // Request a screen refresh
}

/**
 * Draw coordinate axes
 * X axis: red, Y axis: green, Z axis: blue
 */
void drawVelocityVector(){  

    float length = sqrt(velocityx*velocityx + velocityy*velocityy + velocityz*velocityz);
    float normX = velocityx / length; // Normalize velocity vector
    float normY = velocityy / length;
    float normZ = velocityz / length;
    
    glBegin(GL_LINES);
    glColor3f(1, 1, 0); // Red color for velocity vector
    glVertex3f(spherex, spherey, spherez); // Start point (sphere position)
    glVertex3f(spherex +normX, spherey +normY, spherez +normZ); // End point (velocity direction)
    glEnd();
    // Draw arrowhead for velocity vector
    glBegin(GL_TRIANGLES);
    glColor3f(1, 0, 0); // Yellow color for arrowhead

    // Base of the arrowhead
    float arrowLength = 0.2f; // Length of the arrowhead
    float arrowWidth = 0.1f;  // Width of the arrowhead
    glVertex3f(spherex + normX, spherey + normY, spherez + normZ); // Tip of the arrowhead
    glVertex3f(spherex + normX - arrowLength * normX + arrowWidth * upx, 
               spherey + normY - arrowLength * normY + arrowWidth * upy, 
               spherez + normZ - arrowLength * normZ + arrowWidth * upz);
    glVertex3f(spherex + normX - arrowLength * normX - arrowWidth * upx, 
               spherey + normY - arrowLength * normY - arrowWidth * upy, 
               spherez + normZ - arrowLength * normZ - arrowWidth * upz);

    glEnd();
}
void sphereMechanics(){
    // Sphere mechanics can be added here
    float gravity = -9.8f,restitution=0.75f,radius=0.3f,timeStep=0.016f*0.4f; // Gravity acceleration
    float length=6.0f, height=3.0f, depth=6.0f; 
    velocityy += gravity * timeStep; // Update vertical velocity
    spherex += velocityx * timeStep; // Update X position
    spherey += velocityy * timeStep; // Update Y position   
    spherez += velocityz * timeStep; // Update Z position
    if(spherey-radius<=-height){
        spherey=-height+radius;
        velocityy*=-restitution;
    }
    if(spherey+radius>=height){
        spherey=height-radius;
        velocityy*=-restitution;
    }
    if(spherex-radius<=-length){
        spherex=-length+radius;
        velocityx*=-1;
    }
    if(spherex+radius>=length){
        spherex=length-radius;
        velocityx*=-1;
    }
    if(spherez-radius<=-depth){
        spherez=-depth+radius;
        velocityz*=-1;
    }
    if(spherez+radius>=depth){
        spherez=depth-radius;
        velocityz*=-1;
    }
    if(abs(velocityy)<0.0001)velocityy=0;
    float displacementX = velocityx * timeStep;
    float displacementZ = velocityz * timeStep;
    rotationz += (displacementX / radius) * (180.0f / PI); // spin for X movement
    rotationx += (displacementZ / radius) * (180.0f / PI); // spin for Z movement
    rotationAxisy=velocityx*upz-velocityz*upx;
    rotationAxisx=velocityz*upy-velocityy*upz;
    rotationAxisz=velocityy*upx-velocityx*upy;
    rotationAngle+=2.5f;;
}
/**F
 * Draw a colored cube centered at the origin
 * Each face has a different color
 */
void drawStripedSphere(float radius, int slices, int stacks) {
    for (int i = 0; i < stacks; i++) {
        float theta1 = i * M_PI / stacks;
        float theta2 = (i + 1) * M_PI / stacks;

        // Are we in the bottom half or top half?
        bool bottomHalf = (i < stacks / 2);

        glBegin(GL_QUAD_STRIP);
        for (int j = 0; j <= slices; j++) {
            float phi = j * 2.0f * M_PI / slices;

            // Alternate stripes, but flip for top vs bottom
            bool isEven = (j % 2 == 0);
            if (bottomHalf) {
                // Bottom: red-green-red...
                glColor3f(isEven ? 1.0f : 0.0f, isEven ? 0.0f : 1.0f, 0.0f);
            } else {
                // Top: green-red-green...
                glColor3f(isEven ? 0.0f : 1.0f, isEven ? 1.0f : 0.0f, 0.0f);
            }

            float x1 = radius * sinf(theta1) * cosf(phi);
            float y1 = radius * cosf(theta1);
            float z1 = radius * sinf(theta1) * sinf(phi);

            float x2 = radius * sinf(theta2) * cosf(phi);
            float y2 = radius * cosf(theta2);
            float z2 = radius * sinf(theta2) * sinf(phi);

            glVertex3f(x1, y1, z1);
            glVertex3f(x2, y2, z2);
        }
        glEnd();
    }
}

void calculateVelocity(){
    float x = (float)rand() / RAND_MAX;
    float y = (float)rand() / RAND_MAX;
    float z = (float)rand() / RAND_MAX;

    float sum = sqrt(x*x + y*y + z*z);
    x /= sum;
    y /= sum;
    z /= sum;

    velocityx = x * velocity;
    velocityy = y * velocity;
    velocityz = z * velocity;
}
void spherePosition(){
    
    float startX = -6.0f;
    float startZ = -6.0f;
    float tileSize = 12.0f / 20;
    int safeCols = 20 - 2 * (0.3f / tileSize);
    int safeRows = 20 - 2 * (0.3f / tileSize);

    // Pick tile indices avoiding edges
    int newCol = rand() % safeCols + (0.3f / tileSize);
    int newRow = rand() % safeRows + (0.3f / tileSize);
    
    spherex = startX + newCol * tileSize + tileSize / 2.0f;
    spherez = startZ + newRow * tileSize + tileSize / 2.0f;
    spherey = -2.7f; // Sit just above the floor

    //calculateVelocity(); // Calculate initial velocity

}
void drawCheckerboard() {
    int rows = 20, cols = 20;           // 50x50 tiles for nice detail
    float startX = -6.0f;
    float startZ = -6.0f;
    float tileSize = 12.0f / cols;      // Total width is 10 units
    
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            if ((i + j) % 2 == 0)
                glColor3f(1.0f, 1.0f, 1.0f); // white
            else
                glColor3f(0.0f, 0.0f, 0.0f); // black
    
            float x = startX + j * tileSize;
            float z = startZ + i * tileSize;
    
            glBegin(GL_QUADS);
            glVertex3f(x, -3.0f, z);
            glVertex3f(x + tileSize, -3.0f, z);
            glVertex3f(x + tileSize, -3.0f, z + tileSize);
            glVertex3f(x, -3.0f, z + tileSize);
            glEnd();
        }
    }
    
    glPushMatrix();
    glTranslatef(spherex,spherey,spherez); // Position sphere so it sits above the checkerboard
    //glRotatef(rotationx, 1.0f, 0.0f, 0.0f); // Rotate around X axis
    //glRotatef(rotationz, 0.0f, 0.0f, 1.0f); // Rotate around Z axis
    glRotatef(rotationAngle, rotationAxisx, rotationAxisy, rotationAxisz); // Rotate around arbitrary axis
    drawStripedSphere(0.3f, 30, 30);
    glPopMatrix();
    
}
void drawCube()
{
    float width = 6.0f, height = 3.0f, depth = 6.0f; // Cube dimensions
    glBegin(GL_QUADS);
    glColor3f(0.675f,0.455f,0.910f);
    glVertex3f( width,  height, -depth);
    glVertex3f(-width,  height, -depth);
    glVertex3f(-width,  height,  depth);
    glVertex3f( width,  height,  depth);

    // Bottom face (y = -height)
    // glColor3f(1.0f, 0.5f, 0.0f); // Orange
    // glVertex3f( width, -height,  depth);
    // glVertex3f(-width, -height,  depth);
    // glVertex3f(-width, -height, -depth);
    // glVertex3f( width, -height, -depth);

    // Front face (z = +depth)
    glColor3f(0.392f,0.949f,0.949f);
    glVertex3f( width,  height,  depth);
    glVertex3f(-width,  height,  depth);
    glVertex3f(-width, -height,  depth);
    glVertex3f( width, -height,  depth);

    // Back face (z = -depth)
    glColor3f(0.867f,0.922f,0.122f);
    glVertex3f( width, -height, -depth);
    glVertex3f(-width, -height, -depth);
    glVertex3f(-width,  height, -depth);
    glVertex3f( width,  height, -depth);

    // Left face (x = -width)
    glColor3f(0.941f,0.322f,0.322f);
    glVertex3f(-width,  height,  depth);
    glVertex3f(-width,  height, -depth);
    glVertex3f(-width, -height, -depth);
    glVertex3f(-width, -height,  depth);

    // Right face (x = +width)
    glColor3f(0.408f,0.859f,0.298f);
    glVertex3f( width,  height, -depth);
    glVertex3f( width,  height,  depth);
    glVertex3f( width, -height,  depth);
    glVertex3f( width, -height, -depth);

    glEnd();
    drawCheckerboard();
}


/**
 * Draw a pyramid with color gradients
 * Base at y=-1, apex at y=1
 */
void timerFunc(int value)
{
    bladeang+=bladespeed;
    glutPostRedisplay();
    glutTimerFunc(16,timerFunc,0);
}

/**
 * Main function: Program entry point
 */
int main(int argc, char **argv)
{
    // Initialize GLUT
    glutInit(&argc, argv);

    // Configure display mode and window
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(640, 640);
    glutInitWindowPosition(50, 50);
    glutCreateWindow("OpenGL 3D Drawing");

    // Register callback functions
    glutDisplayFunc(display);
    glutReshapeFunc(reshapeListener);
    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);

    spherePosition(); // Initialize sphere position

    // Initialize OpenGL settings
    initGL();

    glutTimerFunc(0, timerFunc, 0); // Start the timer for animation

    // Enter the GLUT event loop
    glutMainLoop();

    return 0;
}