#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>

#define PI 3.1415926
typedef struct Threat {
  double poo_x;
  double poo_y;
  double shooting_height;
  double  corrected_drag_coefficient;
  double muzzleVelocity;
  double ball_m;
  double diameter;
  double theta1;
  double theta2;

}Threat;
Threat threat;

double Atmosphere_drag=0;
typedef struct Bullet {
  char name[50];
  double dragCoefficient;
 
} Bullet;

Bullet bullets[] = {
  {"M43", 0.51,},
  {"SS109", 0.51,},
  {"M855", 0.51, },
  {"MK 262", 0.24,},
  {"Berger Bullets VLD", 0.24,},
  {"Lapua Scenar", 0.24, },
  {"grenade", 0.47, }
};
double flag=0;
double yrut[3] = {0, 0, 0};
double wind[3] = {0, 0, 0};
double air_density=0;
double count=0;

int getRandomBullet() {
  // הגדרת מספר קליעים
  int numBullets = sizeof(bullets) / sizeof(Bullet);
  // בחירת אינדקס אקראי
  int randomIndex = rand() % numBullets;
  // החזרת קליע אקראי
  return randomIndex;
}

double calculateMuzzleVelocity(Bullet bullet) {
  if (bullet.dragCoefficient != 0.47) {
    return (double) rand() / (RAND_MAX / (970 - 350)) + 350;
  } else {
    //רימון
    return (double) rand() / (RAND_MAX / (70 - 20)) + 20;
  }
}

double calculateMasa(Bullet bullet) {
  if (bullet.dragCoefficient != 0.47) {
    // מספר אקראי ממשי בטווח של 0.02 עד 0.110
     return 0.02 + (double) rand() / (double) RAND_MAX * 0.09;
  } else {
    //מספר אקראי ממשי בטווח של 0.3 עד 1.2.
    return 0.3 + (double) rand() / (RAND_MAX + 1.0) * 0.9;
  }
}
double calculate_diameter(Bullet bullet) {
  if (bullet.dragCoefficient != 0.47) {
    //מ1.5 סנטים עד 0.5 סנטים
    return (double) (rand() / (RAND_MAX / (1.5 - 0.5)) + 0.5)/100;
  } else {
    //מגריל ביו 10 סנטים עד סנטים 1
    return (double) (((int) rand() % 10) + 1.0)/100;
  }
}
double Flag_c(){
    return flag;}  
double yrut_x_c(){
    return yrut[0];} 
double yrut_y_c(){
    return yrut[1];} 
double yrut_z_c(){
    return yrut[2];}  
double count_c(){
    return count;}                    
double ball_m_c(){
    return threat.ball_m;}
double corrected_drag_coefficient_c(){
    return threat.corrected_drag_coefficient;}   
double diameter_c(){
    return threat.diameter;}
double muzzleVelocity_c(){
    return threat.muzzleVelocity;}  
double poo_x_c(){
    return threat.poo_x;}
double poo_y_c(){
    return threat.poo_y;} 
double shooting_height_c(){
    return threat.shooting_height;}  
double theta1_c(){
    return threat.theta1;}
double theta2_c(){
    return threat.theta2;} 
 

void *ballistic_trajectory(void* arg){
    clock_t start, end;
    double cpu_time_used;
    start = clock();
    count+=1.0;
    printf("threat number: %0.0lf\n", count);
    double r0[3] = {threat.poo_x, threat.shooting_height, threat.poo_y};  // מיקום ראשוני  (x, y, z)
    double g[3] = {0, 9.81, 0};  
    double m =  threat.ball_m ; 
    double v0 = threat.muzzleVelocity;  
    double point_3[3] = {0, 0, 0};
    int local_count=0;
    yrut[0] = 0;
    yrut[1] = 0;
    yrut[2] = 0;
    //  המרה לרדיאנים
    //זווית אנכית
    double vertical_angle = threat.theta1 * PI / 180;
    //זווית אופקית
    double horizontal_angle = threat.theta2  * PI / 180;
    // חישוב תחילת המומנטום
    double x = m * v0 * cos(horizontal_angle)+wind[0]*m;
    double y = m * v0 * sin(horizontal_angle)+wind[1]*m;
    double z = m * v0 * cos(vertical_angle)+wind[2]*m;
    double p[3] = {x, y, z};  // Momentum (x, y, z)
    // פרמטרים של הסימולצייה
    double rho = air_density;  // Air density
    double A = threat.diameter;
    double t = 0;
    double dt = 0.001;
    // Corrected drag coefficient
    double corrected_drag_coefficient = threat.corrected_drag_coefficient;
    flag=0;
    printf("start c position (x, y, z): (%f, %f, %f)\n", r0[0], r0[1], r0[2]);
    while (abs(r0[0])<150&&r0[1]<150&&r0[1]>-10&&abs(r0[2])<150) {
        //עדכון מהירות
        double v[3] = {p[0] / m, p[1] / m, p[2] / m};  
        local_count++;
        double mag_ball_v = sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2));
        double Fnet[3] = {
            -0.5 * rho * A * corrected_drag_coefficient * ((v[0] / mag_ball_v) * pow(mag_ball_v, 2)) ,
           - (g[1] *m) - 0.5 * rho * A * corrected_drag_coefficient * ((v[1] / mag_ball_v) * pow(mag_ball_v, 2)) ,
            -0.5 * rho * A * corrected_drag_coefficient * ((v[2] / mag_ball_v) * pow(mag_ball_v, 2)) 
        };
           // מיירטים דווקא בנקודה 3 גם כדי שיהיה זמן לחישוב המסלול מהנקודה ה1 לנקודה השלישית בזמן אמת
           //  סופי: אם עובר ב3 מטר אורך רוחב וגובה ומסכן החייל- אז ניירט בנקודה 3
           // ה3 מטר אורך .. חושב כך שזה המרחק שלוקח במקסימום לחייל לרוץ בזמן ריצת הלולאה בפעם 1
            // בהנחה שאנו מחשבים רק את ה4 מטרים אחרי סביבנו כלומר לא יתכן מצב שהמיקום יעבור דרכי לפני שחושבה הנקודה השלישית של הירוט כלומר שהנקודה שלי תהיה לפני השלישית
          if (abs(r0[0])<=3.0&& abs(r0[1]) <=3.0&& abs(r0[2])<=3.0||flag ){
            if (local_count > 3.0) 
            { yrut[0]=point_3[0];
              yrut[1]=point_3[1]; 
              yrut[2]=point_3[2]; 
          flag = 1;
          printf("\ninterception !\n");
           printf("interception point (x, y, z): (%f, %f, %f)\n", yrut[0], yrut[1], yrut[2]);
          return 0;
              }  
             else flag = 1;
        } 
        // עדכן מומנטום ומיקום
        for (int i = 0; i < 3; i++) {
            p[i] =p[i]+ Fnet[i] * dt;
            r0[i] =r0[i]+ (p[i] * dt) / m;
        }
         if(local_count==3.0)
            {
              point_3[0]=r0[0];
              point_3[1]=r0[1]; 
              point_3[2]=r0[2];
        }
        t += dt;  // Update time 
    }
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("%f Time the route will take in seconds\n", cpu_time_used);
    printf("Final position of the route in c (x, y, z): (%f, %f, %f)\n", r0[0], r0[1], r0[2]);
}

void random_()
{
 clock_t start = clock();
	//כדי שכל פעם יגריל לי מספר אחר
	time_t t = time(NULL);
	srand(t);
	double shooting_height, poo_x, poo_y, istabroot;
	//כדי שיאתר לי כל הזמן איומים בלי סוף לכן לולאה אינסופית
  while (true) {
    threat.poo_x = (rand() % (1500 + 1500 + 1) - 1500)/100;
    threat.poo_y = (rand() % (1500 + 1500 + 1) - 1500)/100;
  // אני לא הופכת לערך מוחלט כי כל החישובים מעלים בשנייה שיוצא תמיד חיובי
  // הגרלת מרחק ניתן מהסנסור חייב להתחיל במרחק מינימלי בהגרלה כי נתון הנקודת מוצא
  // פיתגורס
  double Min_distance = (sqrt(pow(threat.poo_x, 2) + pow(threat.poo_y, 2)));
  // חישוב מרחק
  //15 בשביל הסימולצייה
  double distance = (double)rand() / RAND_MAX * (15 - Min_distance) + Min_distance;
  // חישוב גובה מחבל ביחס לגובה חייל
  // יתר בשנייה כפול גובה בשנייה (נעלם כרגע) על כל זה שורש = למרחק
  double hypotenuse = (double)sqrt((threat.poo_x * threat.poo_x) + (threat.poo_y * threat.poo_y));
  
  // if (hypotenuse > distance) {
      
  // } else
   if (hypotenuse == distance) {
      threat.shooting_height = 0;
  } else if (hypotenuse < distance){
      threat.shooting_height = sqrt(pow(distance, 2) - pow(hypotenuse, 2));
  }
    //אם האיום לא מה4 מטר אורך גובה ורוחב שלי תמשיך בהגרלות
    //חושב 4 מטר כי ה3 יכול להיות ממש נקודת ירוט בראש החייל
  if (!(abs(threat.shooting_height)<=7.0&& abs(threat.poo_x) <=7.0&& abs(threat.poo_y)<=7.0)){
        int num_threat=getRandomBullet();
  if(num_threat==6){
         printf("\ngrenade\n");
        }
  else{
       printf("\nbullet\n");
     }
  threat.corrected_drag_coefficient = bullets[num_threat].dragCoefficient*Atmosphere_drag;	
  threat.ball_m =calculateMasa(bullets[num_threat]);
 //הגרלת מספר אקראי בין מהירות תעופת רימון נפץ לבין המהירות הגבוה של קליעים
  threat.muzzleVelocity = calculateMuzzleVelocity(bullets[num_threat]);
  threat.theta1 = rand() % 360;
  threat.theta2 = rand() % 360;
  threat.diameter= calculate_diameter(bullets[num_threat]);
  //שערוך מסלול
  pthread_t tid;
  pthread_create(&tid,NULL,&ballistic_trajectory,NULL);  
  pthread_join(tid,NULL);
  clock_t end = clock();
  double seconds = (double)(end - start) / CLOCKS_PER_SEC;
  sleep(5);
}
}}

void atmosphere_correction( double altitude, double barometer, double temperature, double relative_humidity,double windx,double windy,double airdensity) {
  wind[0]=windx;
  wind[1]=windy;
  air_density=airdensity;
  //#המרה מיחידות מיליבר ליחידות אינטץ עופרת
  barometer = barometer * 0.02953;
  printf("Altitude: %f\n", altitude);
  printf("Barometer: %f\n", barometer);
  printf("Temperature: %f\n", temperature);
  printf("Relative humidity: %f\n\n", relative_humidity);
  double fa = 1/((-4e-15)*pow(altitude, 3) + (4e-10)* pow(altitude, 2)+(-3e-5)*altitude+1);
  double ft = (temperature-(-0.0036 * altitude + 59))/ (459.6 +(-0.0036 * altitude + 59));
  double VPw = (4e-6) * pow(temperature, 3) +-0.0004 * pow(temperature, 2) +   0.0234* temperature +  -0.2517;
  double fr = 0.995 * barometer / (barometer - 0.3783 * relative_humidity * VPw);
  double fp = (barometer-29.53)/29.53;
  double cd = (fa * (1 + ft - fp) * fr);
  Atmosphere_drag=cd;
}










