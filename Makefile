g=6 default1 default1 default1 default1 default1 default1
fps=500
t=10
frameskip=5
all: compile

compile:
	javac sail/sim/Simulator.java

gui:
	java sail.sim.Simulator --verbose --fps ${fps} --gui -g ${g} -t ${t} --frameskip ${frameskip}

run:
	java sail.sim.Simulator -g ${g} -t ${t}

verbose:
	java sail.sim.Simulator -g ${g}  -t ${t} --verbose
