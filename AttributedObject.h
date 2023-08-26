///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  AttributedObject.h
//  ------------------------
//  
//  Base code for rendering assignments.
//
//  Minimalist (non-optimised) code for reading and 
//  rendering an object file
//  
//  We will make some hard assumptions about input file
//  quality. We will not check for manifoldness or 
//  normal direction, &c.  And if it doesn't work on 
//  all object files, that's fine.
//
//	Variant on TexturedObject that stores explicit RGB
//	values for each vertex
//  
///////////////////////////////////////////////////

// include guard for AttributedObject
#ifndef _ATTRIBUTED_OBJECT_H
#define _ATTRIBUTED_OBJECT_H

// include the C++ standard libraries we need for the header
#include <vector>
#include <iostream>
#include <array>
#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <GL/gl.h>
#endif

// include the unit with Cartesian 3-vectors
#include "Cartesian3.h"
// the render parameters
#include "RenderParameters.h"

// define a macro for "not used" flag
#define NO_SUCH_ELEMENT -1

// use macros for the "previous" and "next" IDs
#define PREVIOUS_EDGE(x) ((x) % 3) ? ((x) - 1) : ((x) + 2)
#define NEXT_EDGE(x) (((x) % 3) == 2) ? ((x) - 2) : ((x) + 1)


//half edge vertex structure
struct HE_Vertex {

    Cartesian3 pos;                  //Vertex Position
    int FDedge;                      //first directed edge
    bool is_boundary=false;

    HE_Vertex(){
        pos=Cartesian3();
        FDedge=-1;
        bool is_boundary=false;

    }

    HE_Vertex(Cartesian3 _pos,int _FDedge){
        pos=_pos;
        FDedge=_FDedge;
    }

     HE_Vertex(Cartesian3 _pos){
        pos=_pos;
        FDedge=-1;
    }

     HE_Vertex(int _FDedge){
        pos=Cartesian3();
        FDedge=_FDedge;
    }
};

//half edge vertex structure
struct HE_halfEdge {

    int begin;
    int end;
    int next;
    int prev;
    int face;
    bool isBorder;
    int otherHalf;

    HE_halfEdge() {
        begin = -1;
        end = -1;
        next = -1;
        face = -1;
        prev= -1;
        otherHalf = -1;
        isBorder = false;
    }

    HE_halfEdge(int _begin,int _end,int _next,int _prev,int _face){
        begin=_begin;
        end = _end;
        next=_next;
        prev=_prev;
        face=_face;
    }

};

class AttributedObject
    { // class AttributedObject
    public:
    // vector of vertices
    std::vector<Cartesian3> vertices;

    // vector of colours stored as cartesian triples in float
    std::vector<Cartesian3> colours;
    
    // vector of normals
    std::vector<Cartesian3> normals;
    
    // vector of texture coordinates (stored as triple to simplify code)
    std::vector<Cartesian3> textureCoords;

    // vector for the "first" directed edge for each vertex
    std::vector<long> firstDirectedEdge;

    // vector of faces - doubles as the "to" array for edges
    std::vector<unsigned int> faceVertices;

    // vector of the "other half" of the directed edges
    std::vector<long> otherHalf;

    // centre of gravity - computed after reading
    Cartesian3 centreOfGravity;

    // size of object - i.e. radius of circumscribing sphere centred at centre of gravity
    float objectSize;

    //vertices vector of half edge structure
    std::vector<HE_Vertex> verts;

    //edge vector of half edge structure
    std::vector<HE_halfEdge> halfEdges;

    //to compute all the neighbours of a vertex to update the position of the vertex
    std::vector<int> AdjacentVertex(int i);

    // constructor will initialise to safe values
    AttributedObject();
   
    // read routine returns true on success, failure otherwise
    bool ReadObjectStream(std::istream &geometryStream);

    // write routine
    void WriteObjectStream(std::ostream &geometryStream);

    // routine to render
    void Render(RenderParameters *renderParameters);

    //generate half-edge data structure
    void genHalfEdge();

    void genOtherHalf();

    void genFDEdge();

    void genHEDataStructure();

    //generate
    void genVerts();

    //initialize the texture coordinate vector
    void genTextureCoordinate();

    //to get the normal map
    void genNormals();

    //parameterisation
    void updateBoundaryVerts();

    void updateInteriorVerts();

    void parameterisation();

    Cartesian3 AComputeInteriorVertex(int i);

    //use average weight to update the position of the vertex
    void averageWeight();







    }; // class AttributedObject

// end of include guard for AttributedObject
#endif
