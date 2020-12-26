#pragma once
#include<Model3D.h>
#include <QThread>
#include <iostream>

using namespace std;

class ModelLoder : public QThread
{
	Q_OBJECT


public:

	ModelLoder();
	void setFileName(QString);

	bool loadObjFile(QString);
	bool loadMtlFile(QString);
	

	void loadBinocularModel();


	Model3D getloadedModel();

	~ModelLoder();

private:

	QString obj_file_path;
	QString mtl_file_path;

	Model3D model_3d;

protected:

	void run();

};
