import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib.animation import FuncAnimation

# Constants
# g0 : Acceleration due to gravity (m/s^2)
# R : Specific gas constant for dry air (J/(kg·K))
# T0 : Standard temperature at sea level (K)
# L : Temperature lapse rate (K/m)
# P0 : Standard atmospheric pressure at sea level (Pa)
# gamma : Specific heat ratio for dry air

# User Inputs
print("Welcome to Ben's pocket rocket simulator!\n")
print('⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬛⬛⬛⬛⬛⬛⬛')
print('⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬛⬛⬛⬛⬛⬜⬜⬜⬜⬛')
print('⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬛⬛⬛⬛⬜⬜⬜⬜⬜⬜⬜⬛')
print('⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬛⬛⬛⬜⬜⬜⬜⬜⬜⬜⬜⬜⬛⬛')
print('⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬛⬛⬛⬜⬜⬜⬛⬛⬜⬜⬜⬜⬜⬛⬜')
print('⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬛⬛⬜⬜⬜⬜⬛⬛⬛⬛⬜⬜⬜⬛⬛⬜')
print('⬜⬜⬜⬜⬜⬜⬜⬜⬜⬛⬛⬜⬜⬜⬜⬛⬛🟦🟦⬛⬛⬜⬜⬛⬛⬜')
print('⬜⬜⬜⬜⬛⬛⬛⬛⬛⬛⬜⬜⬜⬜⬜⬛⬛🟦🟦⬛⬛⬜⬛⬛⬜⬜')
print('⬜⬜⬜⬛⬛⬜⬜⬛⬛⬜⬜⬜⬜⬜⬜⬜⬛⬛⬛⬛⬜⬜⬛⬛⬜⬜')
print('⬜⬜⬜⬛⬜⬜⬛⬛⬜⬜⬜⬜⬜⬜⬜⬜⬜⬛⬛⬜⬜⬛⬛⬜⬜⬜')
print('⬜⬜⬛⬛⬜⬛⬛⬜⬜⬜⬜⬜⬛⬛⬜⬜⬜⬜⬜⬜⬛⬛⬜⬜⬜⬜')
print('⬜⬜⬛⬛⬜⬛⬛⬜⬜⬜⬜⬛⬛⬛⬛⬜⬜⬜⬜⬜⬛⬛⬜⬜⬜⬜')
print('⬜⬛⬛⬛⬛⬛⬜⬜⬜⬜⬛⬛⬛⬛⬛⬜⬜⬜⬜⬛⬛⬜⬜⬜⬜⬜')
print('⬜⬛⬛⬛⬛⬛⬛⬜⬜⬛⬛⬛⬛⬛⬜⬜⬜⬜⬛⬛⬜⬜⬜⬜⬜⬜')
print('⬜⬜⬜⬜⬜🟧🟨🟨⬛⬛⬛⬛⬛⬜⬜⬜⬜⬛⬛⬜⬜⬜⬜⬜⬜⬜')
print('⬜⬜🟧🟧🟨🟨🟨🟨⬛⬛⬛⬛⬜⬜⬜⬜⬛⬛⬜⬜⬜⬜⬜⬜⬜⬜')
print('⬜🟥🟥🟧🟨🟨🟨🟨⬛⬛⬛⬜⬜⬜⬜⬛⬛⬛⬜⬜⬜⬜⬜⬜⬜⬜')
print('🟧🟥🟥🟥🟧🟨🟨🟨🟨🟨🟨⬜⬜⬛⬛⬜⬛⬛⬜⬜⬜⬜⬜⬜⬜⬜')
print('⬜⬜⬜🟥🟧🟨🟨🟨🟨🟨🟧⬛⬛⬛⬜⬜⬛⬛⬜⬜⬜⬜⬜⬜⬜⬜')
print('⬜⬜🟥🟥🟥🟧🟨🟨🟧🟨🟥⬛⬛⬜⬜⬜⬛⬛⬜⬜⬜⬜⬜⬜⬜⬜')
print('⬜🟥🟥🟥🟥🟧🟧🟥🟧🟨🟥⬛⬛⬛⬛⬛⬛⬜⬜⬜⬜⬜⬜⬜⬜⬜')
print('🟧🟥🟥🟥🟥🟥🟥🟥🟥🟧🟥⬛⬛⬛⬛⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜')
print('🟥⬜⬜⬜🟥🟥🟥🟧🟥🟧⬜⬛⬛⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜')
print('⬜⬜⬜🟥🟥🟥⬜🟥🟥⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜')
print('⬜⬜⬜🟥🟧⬜⬜🟥🟧⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜')

# User Inputs
print("Please input the following values to continue")

area = np.longdouble(input("Cross-sectional area (m²): "))
mass = np.longdouble(input("Rocket mass (kg): "))
impulse = np.longdouble(input("Mass-Specific impule (s): "))
thrust = np.longdouble(input("Thrust duration (s): "))

A = area # Cross-sectional area (m^2)
rocket_mass = mass  # kg
isp = impulse  # Specific Impulse (s)
thrust_duration= thrust # How long will the engine burn for (s)



time_step = 0.01  # s
end_time = 100.0  # s

# Create a time array
time_values = [i * time_step for i in range(int(end_time / time_step))]


# Define ISA (International Standard Atmosphere) model for atmospheric properties
def atmosphere(altitude):

    # P0 * (1 - L * altitude / T0) ** (g0 / (R * L))
    pressure = 101325.0 * (1 - 0.0065 * altitude / 288.15) ** (9.80665 / (287.05 * 0.0065))
    # T0 - L * altitude
    temperature = 288.15 - 0.0065 * altitude
    # pressure / (R * temperature)
    density = pressure / (287.05 * temperature)
    return density, pressure, temperature

# Function to calculate gravitational acceleration as a function of altitude
def gravity(altitude):
    # g0 + (g0 / ((1 + altitude / 6.371e6) ** 2))
    return 9.80665 + (9.80665 / ((1 + altitude / 6.371e6) ** 2))

# Function to calculate speed of sound as a function of altitude
def speed_of_sound(altitude, temperature):
    # math.sqrt(gamma * R * abs(temperature))
    return math.sqrt(1.4 * 287.05 * abs(temperature))

def drag_coefficient(mach):
    return 0.2 + 0.1 * math.exp(-0.12 * mach)

def rocket_acceleration(y, t, thrust_function, drag_coefficient, rocket_mass):
    altitude, velocity = y
    density, pressure, temperature = atmosphere(altitude)

    # Check if thrust duration has passed
    if t >= thrust_duration:
        thrust = 0.0
    else:
        # Thrust
        thrust = thrust_function(t)

    # Calculate Mach number
    mach = velocity / speed_of_sound(altitude, temperature)

    # Calculate drag
    drag = 0.5 * density * velocity**2 * A * drag_coefficient(mach)

    # Gravity
    gravitational_force = rocket_mass * gravity(altitude)

    # Rocket's mass decreases due to propellant consumption
    if thrust > 0.0:
        mass_change = thrust / (isp * gravity(altitude))
        if rocket_mass > mass_change:
            rocket_mass -= mass_change
        else:
            thrust = rocket_mass * isp * gravity(altitude)  # Consume the remaining mass
    else:
        # If thrust is zero, set acceleration to gravity (free fall)
        acceleration = (-gravity(altitude) - drag) / rocket_mass

    # Calculate acceleration based on thrust and gravity
    if thrust > 0.0:
        acceleration = (thrust - drag - gravitational_force) / rocket_mass
    else:
        # If thrust is zero, set acceleration to gravity (free fall)
        acceleration = (-gravity(altitude) - drag)/ rocket_mass

#    if wrench == '1':
#        tune_velocity = float(input("Velocity tuning parameter: "))
#        velocity = velocity * tune_velocity       
#    elif wrench == '2':
#        tune_acceleration = float(input("Acceleration tuning parameter: "))
#        acceleration = acceleration * tune_acceleration
#    elif wrench == '3':
#        tune_velocity = float(input("Velocity tuning parameter: "))
#        tune_acceleration = float(input("Acceleration tuning parameter: "))
#        acceleration = acceleration * tune_acceleration
#        velocity = velocity * tune_velocity        
#    else:
#        print("Invalid choice. Please select a valid thrust profile (1-3).")


    return [velocity, acceleration]


# Function to simulate rocket and plot altitude vs. time
def rocket_simulation(rocket_mass, thrust_profile, drag_coefficient, time_values):
    altitudes = []
    velocities = []

    y0 = [0.0, 0.0]  # Initial conditions (altitude, velocity)
    
    for i, t in enumerate(time_values):
        y = odeint(rocket_acceleration, y0, [t, t + time_step], args=(thrust_profile, drag_coefficient, rocket_mass))
        altitude, velocity = y[-1]  # Get the values at the end of the time step
        altitudes.append(altitude)
        velocities.append(velocity)
        
        # Update the maximum altitude
        update_max_altitude(altitude, t)

        y0 = [altitude, velocity]  # Update initial conditions for the next time step
        
        # Update the rocket's altitude list
        rocket_altitude[i] = altitude

    return altitudes, velocities

# Initialize a variable to keep track of the maximum altitude
max_altitude = 0.0
max_altitude_time = 0.0

# Function to update maximum altitude
def update_max_altitude(altitude, time):
    global max_altitude, max_altitude_time
    if altitude > max_altitude:
        max_altitude = altitude
        max_altitude_time = time

# Create a list to store the rocket's altitude values
rocket_altitude = [0.0] * len(time_values)

# Define thrust profile functions
# Constant Thrust (ideal rocket):
def constant_thrust(t):
    return max_thrust  # Set max_thrust to your desired constant thrust value

# Linear Decrease in Thrust (single stage):
def linear_decrease_thrust(t):
    if t <= thrust_duration:
        return max_thrust - thrust_rate * t
    else:
        return 0.0

# Stepwise Thrust Profile (multiple-stage):
def stepwise_thrust(t):
    if t < stage1_burn_time:
        return max_thrust_stage1
    elif t < stage1_burn_time + stage2_burn_time:
        return max_thrust_stage2
    else:
        return 0.0

# Thrust Pulse Profile (periodic bursts of thrust):
def thrust_pulse(t):
    if t % pulse_period < pulse_duration:
        return max_thrust_pulse
    else:
        return 0.0

# Choose thrust profile:
print("Thrust Profile Options:")
print("(1) Ideal rocket")
print("(2) Single Stage")
print("(3) Multiple stages")
print("(4) Periodic Burst")
profile = input("Choose your thrust profile: ")

if profile == '1':
    # Constant Thrust (ideal rocket)
    max_thrust = np.longdouble(input("Enter the constant thrust value (N): "))
    thrust_profile = constant_thrust

elif profile == '2':
    # Single Stage (Linear Decrease in Thrust)
    max_thrust = np.longdouble(input("Enter the initial thrust value (N): "))
    thrust_rate = np.longdouble(input("Enter the thrust rate (N/s): "))
    thrust_profile = linear_decrease_thrust

elif profile == '3':
    # Multiple Stages (Stepwise Thrust Profile)
    max_thrust_stage1 = np.longdouble(input("Enter the thrust value for stage 1 (N): "))
    max_thrust_stage2 = np.longdouble(input("Enter the thrust value for stage 2 (N): "))
    stage1_burn_time = np.longdouble(input("Enter the burn time for stage 1 (s): "))
    stage2_burn_time = np.longdouble(input("Enter the burn time for stage 2 (s): "))
    thrust_profile = stepwise_thrust

elif profile == '4':
    # Periodic Burst (Thrust Pulse Profile)
    max_thrust_pulse = np.longdouble(input("Enter the thrust value during pulses (N): "))
    pulse_duration = np.longdouble(input("Enter the duration of each thrust pulse (s): "))
    pulse_period = np.longdouble(input("Enter the time between thrust pulses (s): "))
    thrust_profile = thrust_pulse

else:
    print("Invalid choice. Please select a valid thrust profile (1-4).")


# Continue with the simulation using the selected thrust profile and user-defined parameters
predicted_altitude = rocket_simulation(rocket_mass, thrust_profile, drag_coefficient, time_values)

# Print the maximum altitude
print(f"Apogee: {max_altitude:.2f} meters at time {max_altitude_time:.2f} seconds")

# Create a figure and axis for the animation
fig, ax = plt.subplots(figsize=(8, 6))

# Create the initial plot
altitudes, velocities = rocket_simulation(rocket_mass, thrust_profile, drag_coefficient, time_values)
initial_plot, = ax.plot(time_values[:1], altitudes[:1])

# Calculate accelerations separately
accelerations = [rocket_acceleration([altitudes[i], velocities[i]], t, thrust_profile, drag_coefficient, rocket_mass)[1]
                 for i, t in enumerate(time_values)]

# Variable to keep track of animation state
anim_running = True

#tune = input("Use Tuning Parameters? (y or press any key to continue):")
# Tuning parameters
#if tune == 'y':
#    wrench = input("Tune:\n(1) velocity profile\n(2) acceleration profile \n(3) both\n")
#else:
#    wrench = ''

readout = input("Do you wish to have a numerical readout? (y or press any key to continue): ")

if readout == 'y':
    read = True
else:
    read = False

# Function to toggle animation state
def toggle_animation(event):
    global anim_running
    if event.key == ' ': 
        if anim_running:
            ani.event_source.stop()
        else:
            ani.event_source.start()
        anim_running = not anim_running

# Bind the spacebar key to toggle the animation
fig.canvas.mpl_connect('key_press_event', toggle_animation)


# Create a function to update the animation
def update(frame):
    plt.clf()

    # Main Plot: Altitude vs. Time
    plt.subplot(2, 2, (1, 2))  # Create a subplot for altitude
    plt.plot(time_values[:frame], altitudes[:frame], label='Altitude')
    plt.title('Rocket Altitude vs. Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Altitude (m)')
    plt.grid(True)
    plt.legend()  # Add a legend to this subplot

    # Subplot 1: Velocity vs. Time
    plt.subplot(2, 2, 3)  # Create a subplot for velocity
    plt.plot(time_values[:frame], velocities[:frame], color='orange', label='Velocity')
    plt.title('Rocket Velocity vs. Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Velocity (m/s)')
    plt.grid(True)
    plt.legend()  # Add a legend to this subplot

    # Subplot 2: Acceleration vs. Time
    plt.subplot(2, 2, 4)  # Create a subplot for acceleration
    plt.plot(time_values[:frame], accelerations[:frame], color='green', label='Acceleration')
    plt.title('Rocket Acceleration vs. Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Acceleration (m/s^2)')
    plt.grid(True)
    plt.legend()  # Add a legend to this subplot

    if read == True:
        print(f"Alt: {altitudes[frame]}m, Time: {time_values[frame]}s, A: {accelerations[frame]}, V: {velocities[frame]}")


# Create the animation
ani = FuncAnimation(fig, update, frames=len(time_values), repeat=False, interval=3)

# Display the animation
plt.tight_layout()  # Ensure proper layout of subplots
plt.show()
