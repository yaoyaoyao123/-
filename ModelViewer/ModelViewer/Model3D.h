#pragma once
#pragma execution_character_set("utf-8")
#include <QString>
#include<QVector>
#include<QVector3D>
#include<QOpenGLBuffer>
#include<QOpenGLTexture>
#include<QMatrix4x4>
#define PI 3.1415926

//材质
struct Material
{
	QString mtlname;

	QVector3D Ka;
	QVector3D Kd;
	QVector3D Ks;
	QVector3D Tf;

	float Ni = 0;
	float d = 1;

	//贴图或材质球
	bool value_or_diffmap;
	bool value_or_specmap;

	//贴图路径
	QString map_kd_path;
	QString map_ks_path;

};

//网格
class Mesh
{

public:

	Mesh() = default;

	//顶点
	QVector<double> vps, tps, nps;
	QVector<double> allps;

	//面索引
	QVector<std::tuple<LONGLONG, LONGLONG, LONGLONG>>faceindexs;

	//材质
	Material mtl;

	//缓存对象
	QOpenGLBuffer vbo;
	QOpenGLTexture *map_kd_texture;
	QOpenGLTexture *map_ks_texture;

};

//模型
class Model3D
{

public:

	Model3D();

	//mesh数组
	QVector<Mesh>meshes;

	//材质库
	QVector<Material>mtls;

	//模型大小
	uint mesh_count = 0;
	uint mtl_count = 0;

	//模型变换矩阵
	QMatrix4x4 transform_mat;

	//变换
	void setPitch(double);
	void setRoll(double);
	void setYaw(double);
	void rotate();

	void scale(double);

	void translate();

	void restoreAllTransform();

	~Model3D();

private:

	double scale_factor;

	int pitch;
	int yaw;
	int roll;

	int x_step;
	int y_step;
	int z_step;

};

