#include <maya/MFnPlugin.h> 

//#include "audioNodefmod.h"

#include "audioNode.h"
#include "eqNurbsNode.h"


MStatus initializePlugin( MObject obj )
{ 
	MStatus   status;
	MFnPlugin plugin( obj, "Foliativ", "0.9", "Any");

	status = plugin.registerNode("audioNode", audioNode::id, audioNode::creator, audioNode::initialize);
	status = plugin.registerNode("eqNurbsNode", eqNurbsNode::id, eqNurbsNode::creator, eqNurbsNode::initialize);
	
	return status;
}

MStatus uninitializePlugin( MObject obj )
{
	MStatus   status;
	MFnPlugin plugin( obj );

	status = plugin.deregisterNode(audioNode::id);
	status = plugin.deregisterNode(eqNurbsNode::id);

	
	return status;
}