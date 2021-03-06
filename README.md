# robot_handeye_calibration
Robot handeye calibration is critical for any robot application who has to defined a tool frame on the robot body but don't know the exact transformation matrix from the tool frame to some known frame on the robot.  For instance, the scenario of my concern is to evaluate the robot accuracy using Vicon markers.  The marker plate is kind of arbitrarily attached on the robot end effector.  Therefore, the real transformation needs to be estimated by demonstrating several poses.  This problem has been defined as a general solution of the AX=XB problem.  Note the poses for rotation and translation estimation need to contain enough variations so that an unique solution is achievable.
In this package, 30 samples of the Vicon and Robot movements are save in csv files as an example.
hand_eye_calibration_demo.m provide a simulated example to demonstrate the algorithm.
hand_eye_calibration.m is the main function for homogenous matrix estimation.
