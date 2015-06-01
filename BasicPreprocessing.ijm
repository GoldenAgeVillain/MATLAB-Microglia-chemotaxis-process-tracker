//
//	Script to preprocess chemotaxis data
//
// 	Author: 				Renaud Jolivet, 2015
// 	Modified:				Renaud Jolivet, 2015
// 	Last modification:		2015-02-16
//
// 	Copied and modified from MotilityProcessing.ijm 
//

run("Despeckle", "stack");
run("Z Project...", "projection=[Max Intensity] all");
run("StackReg ", "transformation=[Rigid Body]");
