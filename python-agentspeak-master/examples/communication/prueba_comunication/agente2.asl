/* Create a agente2 agent that will be used to send the message to the agent agente1 */

!start.

+!start : true <- 
    .print("Hola Agente 1, soy Agente 2");
    .send(agente1, achieve, message).

+! message : true <- 
    .print("Adios").