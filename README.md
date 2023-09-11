## Ben's Pocket Rocket Simulator: A Cosmic Adventure 🚀

Welcome to Ben's pocket rocket simulator! This Python script will take you on an exciting journey through the skies, allowing you to explore the possibilities of rocketry. So, let's dive in and explore the world of rocket simulation!

### Requirements 🛠️
Before we get started, please ensure you have the following requirements in place:

Python 3.x
* NumPy
* SciPy
* Matplotlib

If these aren't already installed, you can quickly grab them using the following command:
```{}
pip install numpy scipy matplotlib
```

### Launch Sequence 🚀

Begin your adventure by cloning this repository or downloading the Python script (*rocket_simulator.py*) to your computer.
Open a terminal or command prompt, and navigate to the directory where you've stored the script.
Start the simulation by executing the script with the following command:
```{}
python rocket_simulator.py
```

You'll be prompted to enter various parameters for your rocket simulation, such as:

1. **Cross-sectional area (m²):** The area of your rocket's cross-section in the atmosphere.
2. **Rocket mass (kg):** The initial mass of your rocket.
3. **Mass-Specific impulse (s):** A measure of your rocket engine's efficiency.
4. **Thrust duration (s):** How long your rocket's engine will burn.

## **Main Engine Start** Thrust Profiles: Choose Your Rocket's Personality! 🚀

Select from a range of thrust profiles to match your rocket's unique character:

> **Ideal Rocket (Constant Thrust):** Experience the reliability of an ideal rocket that maintains a consistent thrust throughout the entire journey ensuring your rocket's performance is as steady as the North Star.
> * Set the stage with a constant thrust value in Newtons (N).
>
>Single Stage (Linear Decrease in Thrust): Embark on a rocket adventure with a single stage, where thrust starts strong and gradually tapers off like the excitement of a roller coaster ride. 
> * Define the initial thrust value (N).
> * Adjust the rate at which thrust decreases (N/s).
> 
> Multiple Stages (Stepwise Thrust Profile): Choreograph your rocket's performance with multiple stages, each marked by its unique thrust value (N) and burn time (s). This profile allows you to script a multi-act rocket show as your rocket evolves through each stage.
>  * Each stage has its unique thrust value (N) and burn time (s).
> 
> Periodic Burst (Thrust Pulse Profile): Transform your rocket into a cosmic DJ, delivering pulsating bursts of thrust like a dance floor sensation.
> * Specify the thrust value during each pulse (N).
> * Control the duration of thrust pulses (s).
> * Define the intervals between each pulse (s).
>   

Pick the profile that best suits your rocket's aspirations!


## The Cosmic Show: Experience the Simulation 🌌
The script offers an animated display of your rocket's altitude, velocity, and acceleration throughout the simulation:
> *Toggle the animation on and off by pressing the spacebar, allowing you to closely follow your rocket's journey.*

## The Outcome: Reach for the Stars! 🌟
After the simulation concludes, you'll discover the maximum altitude your rocket reached and the time it reached its peak.
> *You can also choose to enable a numerical readout during the animation to monitor altitude, time, acceleration, and velocity in real-time.*

### Contributing: Join the Adventure! 🌠
If you encounter any issues with this rocket simulator or have ideas for improvements, please feel free to contribute by creating an issue or a pull request on the GitHub repository. Your input is valued and appreciated!
