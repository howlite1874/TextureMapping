///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  AttributedObject.cpp
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

// include the header file
#include "AttributedObject.h"

// include the C++ standard libraries we want
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>

// include the Cartesian 3- vector class
#include "Cartesian3.h"

#define MAXIMUM_LINE_LENGTH 1024
#define REMAP_TO_UNIT_INTERVAL(x) (0.5 + (0.5*(x)))
#define REMAP_FROM_UNIT_INTERVAL(x) (-1.0 + (2.0*(x)))

#define N_ITERATIONS 100000

// constructor will initialise to safe values
AttributedObject::AttributedObject()
    : centreOfGravity(0.0,0.0,0.0)
{ // AttributedObject()
    // force arrays to size 0
    vertices.resize(0);
    colours.resize(0);
    normals.resize(0);
    textureCoords.resize(0);
    firstDirectedEdge.resize(0);
    faceVertices.resize(0);
    otherHalf.resize(0);
    verts.resize(0);
    halfEdges.resize(0);
} // AttributedObject()

// read routine returns true on success, failure otherwise
bool AttributedObject::ReadObjectStream(std::istream &geometryStream)
{ // ReadObjectStream()

    // create a read buffer
    char readBuffer[MAXIMUM_LINE_LENGTH];

    // the rest of this is a loop reading lines & adding them in appropriate places
    while (true)
    { // not eof
        // character to read
        char firstChar = geometryStream.get();

        //         std::cout << "Read: " << firstChar << std::endl;

        // check for eof() in case we've run out
        if (geometryStream.eof())
            break;

        // otherwise, switch on the character we read
        switch (firstChar)
        { // switch on first character
        case '#':       // comment line
            // read and discard the line
            geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
            break;

        case 'v':       // vertex data of some type
        { // some sort of vertex data
            // retrieve another character
            char secondChar = geometryStream.get();

            // bail if we ran out of file
            if (geometryStream.eof())
                break;

            // now use the second character to choose branch
            switch (secondChar)
            { // switch on second character
            case ' ':       // space - indicates a vertex
            { // vertex read
                Cartesian3 vertex;
                geometryStream >> vertex;
                vertices.push_back(vertex);
                //                         std::cout << "Vertex " << vertex << std::endl;
                break;
            } // vertex read
            case 'c':       // c indicates colour
            { // normal read
                Cartesian3 colour;
                geometryStream >> colour;
                colours.push_back(colour);
                //                         std::cout << "Colour " << colour << std::endl;
                break;
            } // normal read
            case 'n':       // n indicates normal vector
            { // normal read
                Cartesian3 normal;
                geometryStream >> normal;
                normals.push_back(normal);
                //                         std::cout << "Normal " << normal << std::endl;
                break;
            } // normal read
            case 't':       // t indicates texture coords
            { // tex coord
                Cartesian3 texCoord;
                geometryStream >> texCoord;
                textureCoords.push_back(texCoord);
                //                         std::cout << "Tex Coords " << texCoord << std::endl;
                break;
            } // tex coord
            default:
                break;
            } // switch on second character
            break;
        } // some sort of vertex data

        case 'f':       // face data
        { // face
            // make a hard assumption that we have a single triangle per line
            unsigned int vertexID;

            // read in three vertices
            for (unsigned int vertex = 0; vertex < 3; vertex++)
            { // per vertex
                // read a vertex ID
                geometryStream >> vertexID;

                // subtract one and store them (OBJ uses 1-based numbering)
                faceVertices.push_back(vertexID-1);
            } // per vertex
            break;
        } // face

            // default processing: do nothing
        default:
            break;

        } // switch on first character

    } // not eof

    // compute centre of gravity
    // note that very large files may have numerical problems with this
    centreOfGravity = Cartesian3(0.0, 0.0, 0.0);

    // if there are any vertices at all
    if (vertices.size() != 0)
    { // non-empty vertex set
        // sum up all of the vertex positions
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
            centreOfGravity = centreOfGravity + vertices[vertex];

        // and divide through by the number to get the average position
        // also known as the barycentre
        centreOfGravity = centreOfGravity / vertices.size();

        // start with 0 radius
        objectSize = 0.0;

        // now compute the largest distance from the origin to a vertex
        for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
        { // per vertex
            // compute the distance from the barycentre
            float distance = (vertices[vertex] - centreOfGravity).length();

            // now test for maximality
            if (distance > objectSize)
                objectSize = distance;

        } // per vertex
    } // non-empty vertex set

    // 	std::cout << "Centre of Gravity: " << centreOfGravity << std::endl;
    // 	std::cout << "Object Size:       " << objectSize << std::endl;

    // return a success code
    return true;
} // ReadObjectStream()

// write routine
void AttributedObject::WriteObjectStream(std::ostream &geometryStream)
{ // WriteObjectStream()
    geometryStream << "# " << (faceVertices.size()/3) << " triangles" << std::endl;
    geometryStream << std::endl;

    // output the vertex coordinates
    geometryStream << "# " << vertices.size() << " vertices" << std::endl;
    for (unsigned int vertex = 0; vertex < vertices.size(); vertex++)
        geometryStream << "v  " << std::fixed << vertices[vertex] << std::endl;

    // output the vertex colours
    geometryStream << "# " << colours.size() << " vertex colours" << std::endl;
    for (unsigned int vertex = 0; vertex < colours.size(); vertex++)
        geometryStream << "vc " << std::fixed << colours[vertex] << std::endl;

    // output the vertex normals
    geometryStream << "# " << normals.size() << " vertex normals" << std::endl;
    for (unsigned int vertex = 0; vertex < normals.size(); vertex++)
        geometryStream << "vn " << std::fixed << normals[vertex] << std::endl;

    // output the vertex coords
    geometryStream << "# " << textureCoords.size() << " vertex tex coords" << std::endl;
    for (unsigned int vertex = 0; vertex < textureCoords.size(); vertex++)
        geometryStream << "vt " << std::fixed << textureCoords[vertex] << std::endl;

    // and the faces
    for (unsigned int face = 0; face < faceVertices.size(); face+=3)
    { // per face
        geometryStream << "f";

        // loop through # of vertices
        for (unsigned int vertex = 0; vertex < 3; vertex++)
        { // per vertex
            geometryStream << " ";
            geometryStream << faceVertices[face+vertex] + 1;
        } // per vertex
        // end the line
        geometryStream << std::endl;
    } // per face

} // WriteObjectStream()

// routine to render
void AttributedObject::Render(RenderParameters *renderParameters)
{ // Render()
    // make sure that textures are disabled
    glDisable(GL_TEXTURE_2D);

    float scale = renderParameters->zoomScale;
    scale /= objectSize;
    // Scale defaults to the zoom setting
    glTranslatef(-centreOfGravity.x * scale, -centreOfGravity.y * scale, -centreOfGravity.z * scale);

    if (renderParameters->useWireframe)
        glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    else
        glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

    Matrix4 tbnMatrix;
    tbnMatrix.SetIdentity();


    // start rendering
    glBegin(GL_TRIANGLES);

    // loop through the faces: note that they may not be triangles, which complicates life
    for (unsigned int face = 0; face < faceVertices.size(); face+=3)
    { // per face
        if(renderParameters->renderNormalMap){
            //normal of three vertices
            Cartesian3 n = normals[face];

            //texture coordinate
            Cartesian3 tc0=textureCoords[faceVertices[face]];
            Cartesian3 tc1=textureCoords[faceVertices[face+1]];
            Cartesian3 tc2=textureCoords[faceVertices[face+2]];

            //position of three vertices
            Cartesian3 pos1(verts[faceVertices[face]].pos);
            Cartesian3 pos2(verts[faceVertices[face+1]].pos);
            Cartesian3 pos3(verts[faceVertices[face+2]].pos);

            Cartesian3 edge1=pos2-pos1;
            Cartesian3 edge2=pos3-pos1;

            Cartesian3 tangent;
            Cartesian3 bitangent;

            //delta UV coordinates
            Cartesian3 dUV1=tc1-tc0;
            Cartesian3 dUV2=tc2-tc0;

            float f = 1.0f / (dUV1.x * dUV2.y - dUV2.x * dUV1.y);

            tangent.x = f * (dUV2.y * edge1.x - dUV1.y * edge2.x);
            tangent.y = f * (dUV2.y * edge1.y - dUV1.y * edge2.y);
            tangent.z = f * (dUV2.y * edge1.z - dUV1.y * edge2.z);
            tangent=tangent.unit();

            bitangent.x = f * (-dUV2.x * edge1.x + dUV1.x * edge2.x);
            bitangent.y = f * (-dUV2.x * edge1.y + dUV1.x * edge2.y);
            bitangent.z = f * (-dUV2.x * edge1.z + dUV1.x * edge2.z);
            bitangent=bitangent.unit();

            auto normal=bitangent.cross(tangent).unit();
            tangent = (tangent - tangent.dot(normal) * normal).unit();
            Cartesian3 tbn[3] = {tangent,bitangent,normal};

            for(size_t row = 0;row<3;row++){
                for(size_t col = 0;col<3;col++){
                    tbnMatrix[row][col] = tbn[row][col];
                }
            }
            tbnMatrix=tbnMatrix.transpose();

            n = tbnMatrix * n;
            n = normal.unit();
            n.x = REMAP_TO_UNIT_INTERVAL(n.x);
            n.y = REMAP_TO_UNIT_INTERVAL(n.y);
            n.z = REMAP_TO_UNIT_INTERVAL(n.z);
            normals[face] = normal;
            normals[face + 1] = normal;
            normals[face + 2] = normal;

        }
        // now do a loop over three vertices
        for (unsigned int vertex = 0; vertex < 3; vertex++)
        { // per vertex
            // set colour using vertex ID

            if(renderParameters->useTexCoords ){
                glColor3f
                        (
                            textureCoords[faceVertices[face+vertex]].x,
                        textureCoords[faceVertices[face+vertex]].y,
                        textureCoords[faceVertices[face+vertex]].z
                        );
            }
            else if(renderParameters->renderNormalMap){

                if(renderParameters->useNormal){
                    glColor3f
                            (
                                normals[faceVertices[face+vertex]].x*255,
                            normals[faceVertices[face+vertex]].y*255,
                            normals[faceVertices[face+vertex]].z*255
                            );
                }
                else{
                    glColor3f
                            (
                                normals[faceVertices[face+vertex]].x,
                            normals[faceVertices[face+vertex]].y,
                            normals[faceVertices[face+vertex]].z
                            );
                }
            }
            else{
                glColor3f
                        (
                            colours[faceVertices[face+vertex]].x,
                        colours[faceVertices[face+vertex]].y,
                        colours[faceVertices[face+vertex]].z
                        );

            }
            if(renderParameters->renderTexture || renderParameters->renderNormalMap){
                glVertex3f
                        (
                            scale * textureCoords[faceVertices[face+vertex]].x,
                        scale * textureCoords[faceVertices[face+vertex]].y,
                        scale * textureCoords[faceVertices[face+vertex]].z
                        );
            }

            else {
                glVertex3f
                        (
                            scale * vertices[faceVertices[face+vertex]].x,
                        scale * vertices[faceVertices[face+vertex]].y,
                        scale * vertices[faceVertices[face+vertex]].z
                        );
            }


        } // per vertex
    } // per face


    // close off the triangles
    glEnd();


    // revert render mode
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

}

void AttributedObject::genHalfEdge()
{
    for(auto i=0;i<faceVertices.size()/3.0f;i++){
        HE_halfEdge e0, e1, e2;
        e0.begin = faceVertices[i * 3];
        e1.begin = faceVertices[i * 3 + 1];
        e2.begin = faceVertices[i * 3 + 2];

        e0.end = faceVertices[i * 3 + 1];
        e1.end = faceVertices[i * 3 + 2];
        e2.end = faceVertices[i * 3];

        e0.face = i;
        e1.face = i;
        e2.face = i;

        e0.next = i * 3 + 1;
        e1.next = i * 3 + 2;
        e2.next = i * 3;

        e0.prev = i * 3 + 2;
        e1.prev = i * 3;
        e2.prev = i * 3 + 1;

        halfEdges.emplace_back(e0);
        halfEdges.emplace_back(e1);
        halfEdges.emplace_back(e2);
    }
}

void AttributedObject::genOtherHalf(){
    for (unsigned int i = 0; i < halfEdges.size(); i++) {

        unsigned int start1 = halfEdges[i].begin;
        unsigned int end1 = halfEdges[i].end;

        for (unsigned int j = 0; j < halfEdges.size(); j++) {

            unsigned int start2 = halfEdges[j].end;
            unsigned int end2 = halfEdges[j].begin;

            if (start1 == start2 && end2 == end1 &&i!=j) {
                halfEdges[i].otherHalf = j;
                break;
            }
        }
    }

}

void AttributedObject::genFDEdge()
{
    for (unsigned int i = 0; i < verts.size(); i++) {
        unsigned int start = i;
        for (unsigned int j = 0; j < halfEdges.size(); j++)
        {
            if (halfEdges[j].begin == i)
            {
                verts[i].FDedge = j;
                break;
            }
        }
    }
}

//All functions for generating data
void AttributedObject::genHEDataStructure()
{
    genHalfEdge();

    genOtherHalf();

    genVerts();

    genFDEdge();

    genTextureCoordinate();

    genNormals();

}

void AttributedObject::genVerts()
{
    for (size_t i=0;i<vertices.size();i++) {
        verts.emplace_back(vertices[i]);
    }
}

void AttributedObject::genTextureCoordinate()
{
    Cartesian3 tmp;
    for (size_t i=0;i<vertices.size();i++) {
        textureCoords.emplace_back(tmp);
    }
}

void AttributedObject::genNormals()
{
    for (size_t i=0;i<faceVertices.size()/3;i++) {
        Cartesian3 e1=vertices[i*3+2]-vertices[i*3];
        Cartesian3 e2=vertices[i*3+1]-vertices[i*3];
        Cartesian3 normal=e1.cross(e2);
        normals.emplace_back(normal);
    }
}

//to use average weight algorithm
std::vector<int> AttributedObject::AdjacentVertex(int i)
{
    //get all vector of adjacent vertex
    std::vector<int> adjacent;
    adjacent.clear();
    //get all vector of adjacent vertex
    int FDedge = verts[i].FDedge;
    int FDedgeID = FDedge;
    if(FDedge==-1){
        return adjacent;
    }
    do
    {
        adjacent.push_back(halfEdges[FDedgeID].end);
        FDedgeID = halfEdges[FDedgeID].otherHalf;
        if (halfEdges[FDedgeID].end == i)
        {
            FDedgeID = halfEdges[FDedgeID].next;
        }
        else { break; }
    } while (FDedgeID != FDedge);
    return adjacent;
}

void AttributedObject::updateBoundaryVerts()
{
    //find all the boundary edges first
    std::vector<unsigned int> boundaryEdge;
    for(unsigned int i =0; i < halfEdges.size();i++){
        if(halfEdges[i].otherHalf==-1){
            halfEdges[i].isBorder=true;
            boundaryEdge.emplace_back(i);
        }
    }

    //make boundary edges conected
    std::vector<unsigned int> OrderedBoundaryEdge;
    OrderedBoundaryEdge.emplace_back(boundaryEdge[0]);
    //for (size_t i=0;i<boundaryEdge.size();i++) {
    while(OrderedBoundaryEdge.size()<boundaryEdge.size()){
        for (size_t j=0;j<boundaryEdge.size();j++) {
            if(halfEdges[OrderedBoundaryEdge[OrderedBoundaryEdge.size()-1]].end==halfEdges[boundaryEdge[j]].begin){
                if(halfEdges[OrderedBoundaryEdge[OrderedBoundaryEdge.size()-1]].end==halfEdges[boundaryEdge[0]].begin){
                    goto mark;
                }
                else {
                    OrderedBoundaryEdge.emplace_back(boundaryEdge[j]);
                    break;
                }
            }

        }

    }

    //convert to conected boundary vertices
mark:std::vector<unsigned int> boundaryVertIndex;
    for (size_t i=0;i<OrderedBoundaryEdge.size();i++) {
        boundaryVertIndex.emplace_back(halfEdges[OrderedBoundaryEdge[i]].begin);
        //        boundaryVertIndex.emplace_back(halfEdges[OrderedBoundaryEdge[i]].end);
        verts[halfEdges[OrderedBoundaryEdge[i]].begin].is_boundary=true;
    }
    verts[halfEdges[OrderedBoundaryEdge[OrderedBoundaryEdge.size()-1]].end].is_boundary=true;


    //update the postion of boundary vertex
    //square boundary
    float sumAll;
    for (size_t i=0;i<boundaryVertIndex.size()-1;i++) {
        //sumAll=sumAll+(verts[boundaryVertIndex[i+1]].pos-verts[boundaryVertIndex[i]].pos).abs();
        sumAll = sumAll+ Cartesian3(verts[boundaryVertIndex[i+1]].pos.x-verts[boundaryVertIndex[i]].pos.x,
                verts[boundaryVertIndex[i+1]].pos.y-verts[boundaryVertIndex[i]].pos.y,
                verts[boundaryVertIndex[i+1]].pos.z-verts[boundaryVertIndex[i]].pos.z
                ).length();
    }

    for (size_t j=0;j<boundaryVertIndex.size();j++) {
        float sumBoundary=0;
        if(j>=1){
            for (size_t k=0;k<j-1;k++) {
                //Cartesian3 tmp=verts[boundaryVertIndex[k+1]].pos-verts[boundaryVertIndex[k]].pos;
                //sumBoundary=sumBoundary+tmp.abs();
                sumBoundary = sumBoundary+ Cartesian3(verts[boundaryVertIndex[k+1]].pos.x-verts[boundaryVertIndex[k]].pos.x,
                        verts[boundaryVertIndex[k+1]].pos.y-verts[boundaryVertIndex[k]].pos.y,
                        verts[boundaryVertIndex[k+1]].pos.z-verts[boundaryVertIndex[k]].pos.z
                        ).length();
            }
        }


        float div=sumBoundary/sumAll;

        if(div>=0&&div<0.25){
            textureCoords[boundaryVertIndex[j]]=Cartesian3(0.5-4*div,1,0);
        }
        else if(div>=0.25&&div<0.5){
            textureCoords[boundaryVertIndex[j]]=Cartesian3(-0.5,2-4*div,0);

        }
        else if(div>=0.5&&div<0.75){
            textureCoords[boundaryVertIndex[j]]=Cartesian3(4*div-2.5,0,0);

        }
        else if(div>=0.75&&div<1){
            textureCoords[boundaryVertIndex[j]]=Cartesian3(0.5,4*div-3,0);

        }

    }
}



void AttributedObject::updateInteriorVerts()
{
    //1.set all interior points to the centre
    Cartesian3 init(0,0.5,0);
    for (size_t i=0;i<verts.size();i++) {
        if(verts[i].is_boundary==false){
            textureCoords[i]=init;
        }
    }

    //2.iterate to get better outcome
    size_t it=200;
    for (size_t i=0;i<it;i++) {
        //2.1 update the position of inner vertices
        averageWeight();


    }
}

//transfer an 3d image into 2d image
void AttributedObject::parameterisation()
{
    updateBoundaryVerts();

    updateInteriorVerts();
}

//use average weight to compute the position of vertex i
Cartesian3 AttributedObject::AComputeInteriorVertex(int i)
{
    std::vector<int> adjacentVerts=AdjacentVertex(i);
    Cartesian3 sum;
    if(adjacentVerts.size()==0){
        return Cartesian3(0.0f,0.5f,0.0f);
    }
    for (size_t j=0;j<adjacentVerts.size();j++) {
        sum=sum+textureCoords[adjacentVerts[j]];
    }
    return sum/adjacentVerts.size();

}

void AttributedObject::averageWeight()
{
    for (size_t j=0;j<verts.size();j++) {
        //just update the position of vertices which are not boundary
        if(verts[j].is_boundary==false){
            textureCoords[j]=AComputeInteriorVertex(j);
        }
    }
}















