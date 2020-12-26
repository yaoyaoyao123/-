#pragma once
#include"ModelLoder.h"
#include<QOpenGLWidget>
#include<QOpenGLFunctions>
#include<QOpenGLShaderProgram>
#include<QOpenGLBuffer>
#include<QOpenGLVertexArrayObject>
#include<QWheelEvent>
#include<QKeyEvent >

class ModelShader : public QOpenGLWidget,protected QOpenGLFunctions
{

	Q_OBJECT

		  
private:

	ModelLoder model_loader;
	QVector<Model3D> model_vec;
	bool model_load_comple;
	int current_select_model = 0;

	QOpenGLShaderProgram m_program;

	QVector<double>horizon_points;
	QOpenGLBuffer h_vbo;
	QOpenGLVertexArrayObject h_vao;
	
	int widget_w, widget_h;
	double fov;

	QVector3D m_eye, m_target, another_vec;
	double dis_to_target;
	double came_yaw, came_pitch;

	QMatrix4x4 camera_mat, project_mat;

public:

	ModelShader(QWidget *parent);

	//模型变换
	void rotateModel();
	void scaleModel(double);
	void translateModel();
	void restoreModelTransform();

	//改变渲染模型
	void setModel(QString);
	void valueForModel();
	void addModel(QString);
	void deleteModel();

protected:

	void initializeGL() override;
	void resizeGL(int w, int h) override;
	void paintGL() override;

	void mousePressEvent(QMouseEvent *event) override;
	void mouseMoveEvent(QMouseEvent *event) override;
	void wheelEvent(QWheelEvent *event) override;
	void keyPressEvent(QKeyEvent *event)override;

	QPoint m_last_pos;


};
