g=7 g1 g2 g3 g4 g5 g6 default1
fps=50
t=50
frameskip=1
S=10
timelimit = 5000
timestep = 0.015
all: compile

compile:
	javac sail/sim/Simulator.java

gui:
	java sail.sim.Simulator -tl ${timelimit} -dt ${timestep} --fps ${fps} --gui -g ${g} -t ${t} --frameskip ${frameskip}

run:
	java sail.sim.Simulator -g ${g} -t ${t} -tl ${timelimit} -dt ${timestep}

verbose:
	java sail.sim.Simulator -g ${g}  -t ${t} --verbose -tl ${timelimit} -dt ${timestep}
