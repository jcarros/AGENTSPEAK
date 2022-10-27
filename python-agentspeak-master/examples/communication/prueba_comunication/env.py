# Create a python code that create a environtment for the agente1 and agente2 conversation
# The environtment is a simple chat

import agentspeak
import agentspeak.runtime
import agentspeak.stdlib

import os

# Create the environment
env = agentspeak.runtime.Environment()

# Create the agents, we want 1 for each agent
with open(os.path.join(os.path.dirname(__file__), "agente2.asl")) as source:
    agents = env.build_agents(source, 2, agentspeak.stdlib.actions)
    
with open(os.path.join(os.path.dirname(__file__), "agente1.asl")) as source:
    agents.append(env.build_agent(source, agentspeak.stdlib.actions))

    
# Run the environment
if __name__ == "__main__":
    env.run()