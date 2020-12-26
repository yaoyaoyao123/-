#include "ModelShader.h"

ModelShader::ModelShader(QWidget *parent):QOpenGLWidget(parent),m_eye(0,0,2),m_target(0,0,0),another_vec(0,1,0),dis_to_target(2)
{

	connect(&model_loader, &QThread::finished, this, &ModelShader::valueForModel);

}

void ModelShader::initializeGL()
{
	m_program.addShaderFromSourceFile(QOpenGLShader::Vertex, "./vsrc.txt");
	m_program.addShaderFromSourceFile(QOpenGLShader::Fragment, "./fsrc.txt");
	m_program.link();

	initializeOpenGLFunctions();
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	glClearColor(0.12, 0.12, 0.12, 1.0);

	//生成参考面
	for (double h_y = -1.0; h_y < 1.0 + 0.00005; h_y += 0.05)
	{
		if (::fabs(h_y) <= 1e-6)
			continue;
		horizon_points << -1.0 << h_y << 0.0
						<< 1.0 << h_y << 0.0;
	}
	for (double h_x = -1.0; h_x < 1.0 + 0.00005; h_x += 0.05)
	{
		if (::fabs(h_x) <= 1e-6)
			continue;
		horizon_points << h_x << -1.0 << 0.0
						<< h_x << 1.0 << 0.0;
	}

	horizon_points << -1.0 << 0.0 << 0.0
		<< 1.0 << 0.0 << 0.0
		<< 0.0 << -1.0 << 0.0
		<< 0.0 << 1.0 << 0.0;


	h_vao.create();
	h_vao.bind();
	h_vbo.create();
	h_vbo.setUsagePattern(QOpenGLBuffer::StaticDraw);
	h_vbo.bind();
	h_vbo.allocate(horizon_points.constData(),horizon_points.count()*sizeof(double));
	m_program.enableAttributeArray(0);
	m_program.setAttributeBuffer(0,GL_DOUBLE,0,3,0);





}

void ModelShader::paintGL()
{

	m_program.bind();
	
	//变换矩阵
	QMatrix4x4 camera_mat;
	camera_mat.lookAt(m_eye,m_target,another_vec);
	project_mat.setToIdentity();
	project_mat.perspective(fov, GLdouble(widget_w) / widget_h, 0.01, 1000.0);

	m_program.setUniformValue("camMatrix", camera_mat);
	m_program.setUniformValue("proMatrix", project_mat);

	//绘制参考面

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	QMatrix4x4 normal_mat;
	normal_mat.setToIdentity();
	m_program.setUniformValue("modMatrix", normal_mat);

	m_program.setUniformValue("model_or_ab", 1);
	h_vao.bind();
	m_program.setUniformValue("color", 1.0, 1.0, 1.0, 0.3);

	glDrawArrays(GL_LINES,0, horizon_points.count() / 3 - 4);
	m_program.setUniformValue("color", 1.0, 0.0, 0.0, 1.0);
	glDrawArrays(GL_LINES, horizon_points.count() / 3 - 4, 2);
	m_program.setUniformValue("color", 0.0, 1.0, 0.0, 1.0);
	glDrawArrays(GL_LINES, horizon_points.count() / 3 - 2, 2);
	h_vao.release();

	if (model_load_comple == true)
	{

		m_program.setUniformValue("model_or_ab",0);
		m_program.setUniformValue("lightPos", 40.0f, 30.0f, 100.0f);
		m_program.setUniformValue("viewPos",m_eye);

		for (uint i = 0; i < model_vec.count(); ++i)
		{

			m_program.setUniformValue("modMatrix",model_vec[i].transform_mat);

			for (uint j = 0; j < model_vec[i].mesh_count; ++j)
			{

				model_vec[i].meshes[j].vbo.create();
				model_vec[i].meshes[j].vbo.bind();
				model_vec[i].meshes[j].vbo.allocate(model_vec[i].meshes[j].allps.constData(), model_vec[i].meshes[j].allps.size()*sizeof(double));

				if (model_vec[i].mtl_count != 0)
				{
					m_program.setUniformValue("nonemtl",0);

					if (model_vec[i].meshes[j].mtl.value_or_diffmap == true && model_vec[i].meshes[j].mtl.value_or_specmap == true)
					{
						m_program.setUniformValue("value_or_diffmap", 1);
						m_program.setUniformValue("value_or_sepcmap", 1);

						glActiveTexture(GL_TEXTURE0);
						m_program.setUniformValue("diff_tex",0);
					
						model_vec[i].meshes[j].map_kd_texture->bind();


						glActiveTexture(GL_TEXTURE1);
						m_program.setUniformValue("spec_tex", 1);

						model_vec[i].meshes[j].map_ks_texture->bind();

						m_program.enableAttributeArray(0);
						m_program.enableAttributeArray(1);
						m_program.enableAttributeArray(2);

						m_program.setAttributeBuffer(0, GL_DOUBLE, (model_vec[i].meshes[j].tps.count()) * sizeof(double) + (model_vec[i].meshes[j].nps.count()) * sizeof(double),3,0);
						m_program.setAttributeBuffer(1, GL_DOUBLE, 0, 2, 0);
						m_program.setAttributeBuffer(2, GL_DOUBLE, (model_vec[i].meshes[j].tps.count()) * sizeof(double), 3, 0);

						glDrawArrays(GL_TRIANGLES, 0, model_vec[i].meshes[j].vps.count() / 3);

						m_program.disableAttributeArray(0);
						m_program.disableAttributeArray(1);
						m_program.disableAttributeArray(2);

						model_vec[i].meshes[j].vbo.release();
						model_vec[i].meshes[j].map_kd_texture->release();
						model_vec[i].meshes[j].map_ks_texture->release();

					}

					else if (model_vec[i].meshes[j].mtl.value_or_diffmap == true && model_vec[i].meshes[j].mtl.value_or_specmap == false)
					{

						m_program.setUniformValue("value_or_diffmap", 1);
						m_program.setUniformValue("value_or_sepcmap", 0);

						glActiveTexture(GL_TEXTURE0);
						m_program.setUniformValue("diff_tex", 0);

						model_vec[i].meshes[j].map_kd_texture->bind();

						m_program.enableAttributeArray(0);
						m_program.enableAttributeArray(1);
						m_program.enableAttributeArray(2);

						m_program.setAttributeBuffer(0, GL_DOUBLE, (model_vec[i].meshes[j].tps.count()) * sizeof(double) + (model_vec[i].meshes[j].nps.count()) * sizeof(double), 3, 0);
						m_program.setAttributeBuffer(1, GL_DOUBLE, 0, 2, 0);
						m_program.setAttributeBuffer(2, GL_DOUBLE, (model_vec[i].meshes[j].tps.count()) * sizeof(double), 3, 0);

						glDrawArrays(GL_TRIANGLES, 0, model_vec[i].meshes[j].vps.count() / 3);

						m_program.disableAttributeArray(0);
						m_program.disableAttributeArray(1);
						m_program.disableAttributeArray(2);

						model_vec[i].meshes[j].vbo.release();
						model_vec[i].meshes[j].map_kd_texture->release();
					}

					else if (model_vec[i].meshes[j].mtl.value_or_diffmap == false && model_vec[i].meshes[j].mtl.value_or_specmap == true)
					{
						m_program.setUniformValue("value_or_diffmap", 0);
						m_program.setUniformValue("value_or_sepcmap", 1);

						glActiveTexture(GL_TEXTURE1);
						m_program.setUniformValue("spec_tex", 1);

						model_vec[i].meshes[j].map_ks_texture->bind();

						m_program.enableAttributeArray(0);
						m_program.enableAttributeArray(1);
						m_program.enableAttributeArray(2);

						m_program.setAttributeBuffer(0, GL_DOUBLE, (model_vec[i].meshes[j].tps.count()) * sizeof(double) + (model_vec[i].meshes[j].nps.count()) * sizeof(double), 3, 0);
						m_program.setAttributeBuffer(1, GL_DOUBLE, 0, 2, 0);
						m_program.setAttributeBuffer(2, GL_DOUBLE, (model_vec[i].meshes[j].tps.count()) * sizeof(double), 3, 0);

						glDrawArrays(GL_TRIANGLES, 0, model_vec[i].meshes[j].vps.count() / 3);

						m_program.disableAttributeArray(0);
						m_program.disableAttributeArray(1);
						m_program.disableAttributeArray(2);

						model_vec[i].meshes[j].vbo.release();
						model_vec[i].meshes[j].map_ks_texture->release();
					}

					else
					{
						m_program.setUniformValue("value_or_diffmap", 0);
						m_program.setUniformValue("value_or_sepcmap", 0);

						m_program.setUniformValue("Kd", model_vec[i].meshes[j].mtl.Kd);
						m_program.setUniformValue("Ks", model_vec[i].meshes[j].mtl.Ks);


						m_program.enableAttributeArray(0);
						m_program.enableAttributeArray(2);

						m_program.setAttributeBuffer(0, GL_DOUBLE, (model_vec[i].meshes[j].tps.count()) * (sizeof(double)) + (model_vec[i].meshes[j].nps.count()) * (sizeof(double)), 3, 0);
						m_program.setAttributeBuffer(2, GL_DOUBLE, (model_vec[i].meshes[j].tps.count()) * (sizeof(double)), 3, 0);

						glDrawArrays(GL_TRIANGLES, 0, model_vec[i].meshes[j].vps.count() / 3);

						m_program.disableAttributeArray(0);
						m_program.disableAttributeArray(2);

						model_vec[i].meshes[j].vbo.release();
					}

				}

				if (model_vec[i].mtl_count == 0)
				{

					m_program.setUniformValue("nonemtl", 1);

					m_program.enableAttributeArray(0);
					m_program.enableAttributeArray(2);

					m_program.setAttributeBuffer(0, GL_DOUBLE, (model_vec[i].meshes[j].tps.count()) * sizeof(double) + (model_vec[i].meshes[j].nps.count()) * sizeof(double), 3, 0);
					m_program.setAttributeBuffer(2, GL_DOUBLE, (model_vec[i].meshes[j].tps.count()) * sizeof(double), 3, 0);

					glDrawArrays(GL_POINTS, 0, model_vec[i].meshes[j].vps.count() / 3);

					m_program.disableAttributeArray(0);
					m_program.disableAttributeArray(2);

					model_vec[i].meshes[j].vbo.release();


				}

			}
		}

	}

	glDisable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);

}

void ModelShader::resizeGL(int w, int h)
{

	widget_w = w;
	widget_h = h;

	fov = 45;
	
}

void ModelShader::mousePressEvent(QMouseEvent *event)
{
	m_last_pos = event->pos();
}

void ModelShader::mouseMoveEvent(QMouseEvent *event)
{

	int dx = event->x() - m_last_pos.x();
	int dy = event->y() - m_last_pos.y();

	if (event->buttons() & Qt::LeftButton)
	{
		static uint mousemovecount = 0;
		mousemovecount++;
		if (mousemovecount == 1)
		{
			dx = 0;
			dy = 0;
		}

		came_pitch = came_pitch + 0.04*dy;

		if (came_pitch >= 89.9f)
			came_pitch = 89.9f;

		if (came_pitch <= -89.9f)
			came_pitch = -89.9f;

		m_eye.setY(dis_to_target *cos(came_pitch *PI / 180));
		m_eye.setZ(dis_to_target *sin(came_pitch *PI / 180));

		came_yaw = came_yaw + 0.04*dx;
		m_eye.setY(dis_to_target *cos(came_yaw *PI / 180));
		m_eye.setX(dis_to_target *sin(came_yaw *PI / 180));

		m_eye.normalize();
		m_eye = dis_to_target * m_eye;
		another_vec = QVector3D(0, 0, 1);
		m_target = QVector3D(0, 0, 0);

		m_last_pos = event->pos();

		update();
	}

}

void ModelShader::wheelEvent(QWheelEvent *event)
{
	if (event->delta() > 0)
	{
		m_eye = m_eye + 0.03*m_eye;
	}
	else
	{
		m_eye = m_eye - 0.03*m_eye;
	}

	dis_to_target = m_eye.length();
	event->accept();
	update();


}

void ModelShader::keyPressEvent(QKeyEvent *event)
{

	if (event->key() == Qt::Key_Up)
	{
		model_vec[0].transform_mat.translate(0.0, 0.0, 0.1);
	}
	if (event->key() == Qt::Key_Down)
	{
		model_vec[0].transform_mat.translate(0.0, 0.0, -0.1);
	}
	if (event->key() == Qt::Key_Left)
	{
		model_vec[0].transform_mat.translate(0.1, 0.0, 0.0);
	}
	if (event->key() == Qt::Key_Right)
	{
		model_vec[0].transform_mat.translate(-0.1, 0.0, 0.0);
	}
	if (event->key() == Qt::Key_W)
	{
		model_vec[0].transform_mat.translate(0.0, 0.1, 0.0);
	}
	if (event->key() == Qt::Key_E)
	{
		model_vec[0].transform_mat.translate(0.0, -0.1, 0.0);
	}



	if (event->key() == Qt::Key_A)
	{
		model_vec[0].transform_mat.rotate(2.0, 1.0, 0.0, 0.0);
	}
	if (event->key() == Qt::Key_S)
	{
		model_vec[0].transform_mat.rotate(-2.0, 1.0, 0.0, 0.0);
	}
	if (event->key() == Qt::Key_D)
	{
		model_vec[0].transform_mat.rotate(2.0, 0.0, 1.0, 0.0);
	}
	if (event->key() == Qt::Key_F)
	{
		model_vec[0].transform_mat.rotate(-2.0, 0.0, 1.0, 0.0);
	}
	if (event->key() == Qt::Key_G)
	{
		model_vec[0].transform_mat.rotate(2.0, 0.0, 0.0, 1.0);
	}
	if (event->key() == Qt::Key_H)
	{
		model_vec[0].transform_mat.rotate(-2.0, 0.0, 0.0, 1.0);
	}

	if (event->key() == Qt::Key_Q)
	{
		model_vec[0].transform_mat.scale(0.3);
	}


	update();

}

void ModelShader::setModel(QString filename)
{

	model_load_comple = 0;

	model_loader.setFileName(filename);

	model_loader.start();

}



void ModelShader::valueForModel()
{

	model_vec << model_loader.getloadedModel();

	model_load_comple = 1;

	//model_vec[0].meshes[0].map_kd_texture = new QOpenGLTexture(QImage(QString(model_vec[0].meshes[0].mtl.map_kd_path)).mirrored());

	//rotateModel();
	//translateModel();


	update();
	

}



void ModelShader::scaleModel(double scale_factor )
{

	model_vec[0].transform_mat.scale(scale_factor);
	update();

}


void ModelShader::rotateModel()
{

	//model_vec[0].transform_mat.rotate(90,1,0,0);
	model_vec[0].transform_mat.rotate(90.0, 0.0, 1.0, 0.0);
	update();

}

void ModelShader::translateModel()
{
	model_vec[0].transform_mat.translate(0.0,0.0,0.5);
	
}