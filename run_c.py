import vpython
import requests
import geocoder
import threading
import time
import math
import os

from vpython import *
from ctypes import cdll
from ctypes import c_double
libObject= cdll.LoadLibrary('./liberary.so')

#(הסנסור) מחיקת קובץ האיומים בהדלקת התוכנית
if os.path.exists('data.txt'):
    os.remove('data.txt')
#מזג אוויר
atmosphere_correction=libObject.atmosphere_correction
api_key = "f331e7326602c88de4b4f42982ded2f4"
# מביא מיקום נוכחי
g = geocoder.ip('me')
latitude = g.latlng[0]
longitude = g.latlng[1]
# Build the API request URL
url = f"https://api.openweathermap.org/data/2.5/weather?lat={latitude}&lon={longitude}&appid={api_key}&units=metric"
# Send the API request and get the response
response = requests.get(url)
if response.status_code != 200:
   print(f"Error: API key login error: {response.status_code}")
data = response.json()
PI=3.14159
wind_speed = data["wind"]["speed"]
wind_direction = data["wind"]["deg"]
wind_direction_rad = wind_direction*PI/180
wind_x= wind_speed * cos(wind_direction_rad)
wind_y=wind_speed * sin(wind_direction_rad)
wind = vpython.vector(wind_x, wind_y, 0)
print(f"wind{wind}")
windx=c_double(wind_x)
windy=c_double(wind_y)

KELVIN_OFFSET = 273.15  # קבוע להמרת צלזיוס לקלווין
LAPSE_RATE = 0.0065     # שיעור הירידה בטמפרטורה בקלווין למטר
SEA_LEVEL_PRESSURE = 1013.25  # לחץ אטמוספירי סטנדרטי בגובה פני הים בפסקל
EXPONENT = 0.1903       # חזקה בנוסחת לחץ
# נתוני מזג אוויר
temperature_celsius = data["main"]["temp"]
pressure_hPa = data["main"]["pressure"]
temperature_kelvin = temperature_celsius + KELVIN_OFFSET
pressure_ratio = pressure_hPa / SEA_LEVEL_PRESSURE
altitude = c_double(temperature_kelvin / LAPSE_RATE * (1 - pressure_ratio ** EXPONENT))
#גם אם זה לא תוצאה מדוייקת משתמשים בחישוב זה גם במפות וגם משנה לי הגובה בשביל חישובי צפיפות
#התנגדות אוויר ולא מעבר ואם מציינים לי גובה מסויים והוא מתאים לשאר תנאי מזג אוויר =מתאים לחישוב
barometer = c_double(data["main"]["pressure"])
temperature = c_double(data["main"]["temp"])
relative_humidity = c_double(data["main"]["humidity"])
#משתנה לחישוב התנגדות אוויר
R = 287.05  # Ideal gas constant for air in J/(kg·K)
# Convert units
pressure_Pa = data["main"]["pressure"] * 100  # Convert pressure from hPa to Pa
temperature_K = data["main"]["temp"] + KELVIN_OFFSET # Convert temperature from Celsius to Kelvin
# Calculate virtual temperature Tv (in Kelvin)
Tv = temperature_K / (1 - (0.378 * data["main"]["humidity"] / 100))
# Calculate air density (in kg/m^3)
air_density = pressure_Pa / (R * Tv)
rho=air_density
print(f"rho{rho}")
airdensity=c_double(air_density)
libObject.atmosphere_correction( altitude, barometer, temperature, relative_humidity,windx,windy,airdensity,)
#

poo_x=0
poo_y=0
count=0
shooting_height=0
corrected_drag_coefficient=0
muzzleVelocity=0
check=0
ball_m=0
diameter=0
theta1=0
theta2=0
flag=0
yrut_z=0
yrut_y=0
yrut_x=0 

def func_ball_m():
    global ball_m
    ball_m = libObject.ball_m_c
    ball_m.restype = c_double
    ball_m =  (libObject.ball_m_c())
    print(f"ball_m {ball_m}")

def func_count():
    global count
    count = libObject.count_c
    count.restype = c_double
    count =  (libObject.count_c())
     
def func_flag():
    global flag
    flag = libObject.Flag_c
    flag.restype = c_double
    flag =  (libObject.Flag_c())

def func_poo_x():
    global poo_x
    poo_x = libObject.poo_x_c
    poo_x.restype = c_double
    poo_x =libObject.poo_x_c()
    print(f"poo_x {poo_x}")

def func_poo_y():
    global poo_y
    poo_y = libObject.poo_y_c
    poo_y.restype = c_double
    poo_y = libObject.poo_y_c()
   

def func_yrut_x():
    global yrut_x
    yrut_x = libObject.yrut_x_c
    yrut_x.restype = c_double
    yrut_x = libObject.yrut_x_c()


def func_yrut_y():
    global yrut_y
    yrut_y = libObject.yrut_y_c
    yrut_y.restype = c_double
    yrut_y =  (libObject.yrut_y_c())
  

def func_yrut_z():
    global yrut_z
    yrut_z = libObject.yrut_z_c
    yrut_z.restype = c_double
    yrut_z =  (libObject.yrut_z_c())


def func_shooting_height():
    global shooting_height
    shooting_height = libObject.shooting_height_c
    shooting_height.restype = c_double
    shooting_height = libObject.shooting_height_c()
    print(f"shooting_height {shooting_height}")

def func_corrected_drag_coefficient():
    global corrected_drag_coefficient
    corrected_drag_coefficient = libObject.corrected_drag_coefficient_c
    corrected_drag_coefficient.restype = c_double
    corrected_drag_coefficient =  (libObject.corrected_drag_coefficient_c())
    print(f"corrected_drag_coefficient {corrected_drag_coefficient}")
   
def func_muzzleVelocity():
    global muzzleVelocity
    muzzleVelocity = libObject.muzzleVelocity_c
    muzzleVelocity.restype = c_double
    muzzleVelocity =  (libObject.muzzleVelocity_c())
    print("muzz=",muzzleVelocity)
    
def func_diameter():
    global diameter
    diameter = libObject.diameter_c
    diameter.restype = c_double
    diameter =  (libObject.diameter_c())
    print(f"diameter {diameter}")

def func_theta1():
    global theta1
    theta1 = libObject.theta1_c
    theta1.restype = c_double
    theta1 =  (libObject.theta1_c())
    print(f"horizontal_angle {theta1}")

def func_theta2():
    global theta2
    theta2 = libObject.theta2_c
    theta2.restype = c_double
    theta2 =  (libObject.theta2_c())
    print(f"vertical_angle {theta2}")

def simulation():
   # start_time = time.time()  # מקבל את הזמן הנוכחי לפני הרצת התוכנית
    loc_corrected_drag_coefficient=corrected_drag_coefficient
    light_beam = vpython.cylinder(pos=vector(0,0,0), axis=vector(0,0,0), radius=0.3,  color=color.green)
    ball = vpython.sphere(pos=vector(poo_x, shooting_height, poo_y), radius=.9, color=color.red, make_trail=True)
    r0 = ball.pos
    local_rho=rho
    local_ball_m=ball_m
    local_muzzleVelocity=muzzleVelocity
    print("position start in simulation ", r0)
    PI = 3.1415926
    #זווית אנכית
    vertical_angle=theta1*PI/180
    #זווית אופקית
    horizontal_angle = theta2 * PI / 180
    x = math.cos(horizontal_angle) * local_ball_m * local_muzzleVelocity
    y = math.sin(horizontal_angle) *local_ball_m* local_muzzleVelocity
    z= math.cos(vertical_angle) * local_ball_m * local_muzzleVelocity
    #מומנטום =כח תנע
    ball.p=vpython.vector(x,y,z)+wind*local_ball_m
    local_diameter=diameter
    t = 0
    dt = 0.001
    gravity = 9.81  
    local_flag=flag
    local_yrut_z=yrut_z
    local_yrut_y=yrut_y
    local_yrut_x=yrut_x
    while abs(ball.pos.x)<150 and ball.pos.y<150 and ball.pos.y>-10 and abs(ball.pos.z)<150:
        rate(1000)
        if local_flag==1 and t==0.003:
            light_beam.axis=vector(local_yrut_x,local_yrut_y,local_yrut_z)
            local_flag=0 
            current_time = time.strftime("%H:%M")
            sleep(5)
            ball.size=vector(0,0,0)
            light_beam.size =vector(0,0,0)
            ball.clear_trail()
            with open('data.txt', 'a') as file:
                file.write(f"{local_yrut_x*100}\n")
                file.write(f"{local_yrut_y*100}\n")
                file.write(f"{local_yrut_z*100}\n")
                file.write(current_time)
                file.write("\n")
            return
        #mag גודל וקטור תנע
        #norm ווקטור היחידה
        ball.v = ball.p / local_ball_m  
        Fnet_x=-0.5 * local_rho * local_diameter * loc_corrected_drag_coefficient * ((ball.v.x/(mag(ball.v)))*(mag(ball.v)**2))#norm* mag
        Fnet_y = -(gravity*local_ball_m)- 0.5 * local_rho * local_diameter * loc_corrected_drag_coefficient * ((ball.v.y/(mag(ball.v)))*(mag(ball.v)**2))
        Fnet_z = -0.5 * local_rho * local_diameter * loc_corrected_drag_coefficient *((ball.v.z/(mag(ball.v)))*(mag(ball.v)**2))
        Fnet = vpython.vector(Fnet_x, Fnet_y, Fnet_z)
        ball.p = ball.p + Fnet * dt
        ball.pos = ball.pos + ((ball.p * dt) / local_ball_m)
        t += dt 
        
   # end_time = time.time()  # מקבל את הזמן הנוכחי לאחר הרצת התוכנית
   # בלא תמיד מהירות מדוייקת ביגלל התנגדות אוויר
   # print(f"time py{end_time - start_time} ") 
    print("position end py= ",ball.pos)
    sleep(3)
    ball.size=vector(0,0,0)
    ball.clear_trail()


scene = canvas(width=1300, height=800)  # גודל מותאם אישית לחלון
scene.fullscreen = True  # מציג את הסימולציה במצב מסך מלא  
scene.background = vector(0.8, 0.9, 1.0) 
soldier_color = vector(0.33, 0.42, 0.18) 
ground_color = vector(0.95, 0.85, 0.7)   
head =sphere(pos=vector(0, 0, 0), radius=3, color=soldier_color)
top_half = cylinder(pos=vector(0, -.5, 0), axis=vector(0, -10, 0), radius=2, color=soldier_color)
ground = box(pos=vector(0,-10, 0), size=vector(300, .4, 300), color=ground_color)
#קבלת איומים מהסנסור
t = threading.Thread(target=libObject.random_, )
t.start() 
#הצגת סימולציית אמת
while(1): 
  func_count()
  if (count!=check):  
    check=count
    t1 = threading.Thread(target=func_ball_m)
    t2 = threading.Thread(target=func_poo_x)
    t3 = threading.Thread(target=func_poo_y)
    t4 = threading.Thread(target=func_shooting_height)
    t5 = threading.Thread(target=func_corrected_drag_coefficient)
    t6 = threading.Thread(target=func_muzzleVelocity)
    t7 = threading.Thread(target=func_diameter)
    t8 = threading.Thread(target=func_theta1)
    t9 = threading.Thread(target=func_theta2)
    t1.start()
    t2.start()
    t3.start()
    t4.start()
    t5.start()
    t6.start()
    t7.start()
    t8.start()
    t9.start()

    t1.join()
    t2.join()
    t3.join()
    t4.join()
    t5.join()
    t6.join() 
    t7.join()
    t8.join()
    t9.join()
    
    f=threading.Thread(target=func_flag)
    f.start()
    xxx=threading.Thread(target=func_yrut_x)
    yyy=threading.Thread(target=func_yrut_y)
    zzz=threading.Thread(target=func_yrut_z)
    xxx.start()
    yyy.start()
    zzz.start()
    route = threading.Thread(target=simulation, )
    route.start()    
