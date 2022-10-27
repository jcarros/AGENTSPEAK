/* Create a agentspeak agent that recibe a message from the agente2 agent and send a message to the agente2 agent */

+! message <- 
    .print("Hola agente2, soy agente1");
    .send(agente2, achieve, message).
