#include <GL/glut.h>
#include <math.h>
#include <time.h>
const float PI =3.14159265358979323846;
float secAngle, minAngle, hourAngle;

void initTime() {
    time_t rawtime = time(NULL);
    struct tm *timeinfo = localtime(&rawtime);

    secAngle = timeinfo->tm_sec * 6.0f;
    minAngle = timeinfo->tm_min * 6.0f + timeinfo->tm_sec * 0.1f;
    hourAngle = (timeinfo->tm_hour % 12) * 30.0f + timeinfo->tm_min * 0.5f;
}

void drawCircle(float radius) {
    glBegin(GL_LINE_LOOP);
    for (int i = 0; i < 360; i++) {
        float theta = i * M_PI / 180.0f;
        glVertex2f(radius * cos(theta), radius * sin(theta));
    }
    glEnd();
}

void drawHourTicks() {
    for (int i = 0; i < 60; ++i) {
        float angleDeg = i * 6.0f; // 360/12
        float angleRad = angleDeg * (PI / 180.0f);

        float x1 = 0.545f * cos(angleRad);
        float y1 = 0.545f * sin(angleRad);
        float temp=0.58f;
        float x2 = temp * cos(angleRad);
        float y2 = temp * sin(angleRad);
        if(i%5==0){
            x1 = 0.48f * cos(angleRad);
            y1 = 0.48f * sin(angleRad);
        }
    
        glLineWidth(2.0f);
        glBegin(GL_LINES);
            glVertex2f(x1, y1);
            glVertex2f(x2, y2);
        glEnd();
    }
}

void drawHand(float angle, float length, float width) {
    glPushMatrix();
    glRotatef(-angle, 0, 0, 1);
    glLineWidth(width);
    glBegin(GL_LINES);
        glVertex2f(0.0f, 0.0f);
        glVertex2f(0.0f, length);
    glEnd();
    glPopMatrix();
    
}

void drawMovingSquare(float angle, float radius, float size,float r, float g, float b) {
    glPushMatrix();
    glRotatef(-angle, 0, 0, 1);      // Rotate clockwise
    glTranslatef(0.0f, radius, 0.0f); // Move to the edge (tip of hand)
    glColor3f(r, g, b);     // Green box
    glBegin(GL_QUADS);
        glVertex2f(-size/2, -size/2);
        glVertex2f(size/2, -size/2);
        glVertex2f(size/2, size/2);
        glVertex2f(-size/2, size/2);
    glEnd();
    glPopMatrix();
}


void display() {
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();

    glColor3f(1, 1, 1);
    drawHourTicks();

    // Draw clock border
    glColor3f(1, 1, 1);
    drawCircle(0.6f);

    

    // Draw hour hand
    glColor3f(1, 1, 1);
    drawHand(hourAngle, 0.25f, 6.0f);
    drawMovingSquare(hourAngle, 0.6f, 0.03f,1,1,1); // Draw square at the end of hour hand

    // Draw minute hand
    glColor3f(1, 1, 1);
    drawHand(minAngle, 0.33f, 4.0f);
    drawMovingSquare(minAngle, 0.6f, 0.03f,1,1,1); // Draw square at the end of minute hand

    // Draw second hand
    glColor3f(1, 0, 0);
    drawHand(secAngle, 0.4f, 2.0f);
    drawMovingSquare(secAngle, 0.6f, 0.03f,1,0,0); // Draw square at the end of second hand

    glutSwapBuffers();
}

void timer(int value) {
    // Update angles smoothly
    secAngle += 0.1f;
    if (secAngle >= 360.0f) {
        secAngle -= 360.0f;
    }

    minAngle += 0.1f / 60.0f;       // 6° per min = 0.1° per sec
    if (minAngle >= 360.0f) {
        minAngle -= 360.0f;
    }

    hourAngle += 0.1f / 3600.0f;    // 30° per hour = 0.0083° per sec
    if (hourAngle >= 360.0f) {
        hourAngle -= 360.0f;
    }

    glutPostRedisplay();
    glutTimerFunc(16, timer, 0); // ~60 FPS
}

void reshape(int w, int h) {
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-1, 1, -1, 1);
    glMatrixMode(GL_MODELVIEW);
}

int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(500, 500);
    glutCreateWindow("Smooth Analog Clock");

    initTime();

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutTimerFunc(0, timer, 0);

    glClearColor(0, 0, 0, 1); // Black background

    glutMainLoop();
    return 0;
}
