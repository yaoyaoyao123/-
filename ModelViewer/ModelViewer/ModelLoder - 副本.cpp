#include "ModelLoder.h"
#include<QFileInfo>

using namespace cv;

ModelLoder::ModelLoder()
{
}

void ModelLoder::run()
{

	loadBinocularModel();

}

bool ModelLoder::loadObjFile(QString filename)
{

	QVector<double>vertexPoints, texturePoints, normalPoints;

	if (filename.mid(filename.lastIndexOf('.')) != ".obj")
	{
		qDebug() << "file is not a obj file";
		return false;

	}

	QFile objfile(filename);

	if (objfile.isOpen()) objfile.close();

	if (!objfile.open(QIODevice::ReadOnly))
	{
		qDebug() << "open" << filename << "filed";
		return false;
	}
	else
		qDebug() << "open sucess";

	while (!objfile.atEnd())
	{
		QByteArray linedate = objfile.readLine();

		linedate = linedate.trimmed();

		if (linedate == "")
			continue;

		QList<QByteArray> str_values = linedate.split(' ');
		str_values.removeAll("");
		QString data_type = str_values.takeFirst();

		if (data_type == "mtllib")
		{
			QFileInfo mtl(filename);
			mtl_file_path = mtl.path() + "/" + str_values[0];
		}

		else if(data_type == "v")
		{
			for (int i = 0; i < 3; i++)
			{
				if (str_values[i] != "")
					vertexPoints.push_back(str_values[i].toDouble());
			}
		}
		else if (data_type == "vt")
		{
			for (int i = 0; i < 2; i++)
			{
				if (str_values[i] != "")
					texturePoints.push_back(str_values[i].toDouble());
			}
		}
		else if (data_type == "vn")
		{
			for (int i = 0; i < 3; i++)
			{
				if (str_values[i] != "")
					normalPoints.push_back(str_values[i].toDouble());
			}
		}
		else if (data_type == "usemtl")
		{
			model_3d.mesh_count++;
			Mesh meshi;
			model_3d.meshes << meshi;
			model_3d.meshes[model_3d.mesh_count - 1].mtl.mtlname = str_values[0];

		}
		else if (data_type == "f")
		{

			if (str_values.size() % 3 != 0)
			{
				if (str_values.size() == 4)
				{
					str_values.push_back(str_values.at(0));
					str_values.push_back(str_values.at(2));
				}
				else if (str_values.size() == 5)
				{
					str_values.insert(1, str_values.at(4));
					str_values.insert(2, str_values.at(2));
					str_values.insert(3, str_values.at(1));
					str_values.insert(5, str_values.at(5));
				}
				else if (str_values.size() == 7)
				{
					str_values.insert(3, str_values.at(1));
					str_values.insert(4, str_values.at(2));

				}
				else if (str_values.size() == 8)
				{
					str_values.insert(3, str_values.at(1));
					str_values.insert(4, str_values.at(2));
					str_values.insert(9, str_values.at(6));
					str_values.insert(10, str_values.at(8));
				}
				else
				{
					qDebug() << str_values.size() << str_values;
					continue;
				}
			}

			if (model_3d.meshes.count() <= 0)
			{
				return false;
			}

			std::transform(str_values.begin(), str_values.end(), std::back_inserter(model_3d.meshes[model_3d.mesh_count - 1].faceindexs), [](QByteArray &str)
			{

				QList<QByteArray>intStr = str.split('/');
				return std::make_tuple(intStr.at(0).toLongLong(), intStr.at(1).toLongLong(), intStr.at(2).toLongLong());

			});

		}
		else
			continue;

	}

	for (int i = 0; i < model_3d.mesh_count; ++i)
	{

		for (auto &verFaceInfo:model_3d.meshes[i].faceindexs)
		{
			LONGLONG vIndexs = std::get<0>(verFaceInfo);
			LONGLONG tIndexs = std::get<1>(verFaceInfo);
			LONGLONG nIndexs = std::get<2>(verFaceInfo);

			model_3d.meshes[i].vps << vertexPoints.at(vIndexs * 3 - 3);
			model_3d.meshes[i].vps << vertexPoints.at(vIndexs * 3 - 2);
			model_3d.meshes[i].vps << vertexPoints.at(vIndexs * 3 - 1);

			if (texturePoints.size() != 0 && tIndexs !=0)
			{
				model_3d.meshes[i].tps << texturePoints.at(tIndexs * 2 - 2);
				model_3d.meshes[i].tps << texturePoints.at(tIndexs * 2 - 1);
			}

			model_3d.meshes[i].nps << normalPoints.at(nIndexs * 3 - 3);
			model_3d.meshes[i].nps << normalPoints.at(nIndexs * 3 - 2);
			model_3d.meshes[i].nps << normalPoints.at(nIndexs * 3 - 1);

		}

		model_3d.meshes[i].allps << model_3d.meshes[i].tps << model_3d.meshes[i].nps << model_3d.meshes[i].vps;

	}

	vertexPoints.clear();
	texturePoints.clear();
	normalPoints.clear();
	objfile.close();
	return true;

}


bool ModelLoder::loadMtlFile(QString filename)
{
	QVector<double>vertexPoints, texturePoints, normalPoints;

	if (filename.mid(filename.lastIndexOf('.')) != ".mtl")
	{
		qDebug() << "file is not a mtl file";
		return false;
	}

	QFile mtlfile(filename);

	if (mtlfile.isOpen()) mtlfile.close();

	if (!mtlfile.open(QIODevice::ReadOnly))
	{
		qDebug() << "open" << filename << "filed";
		return false;
	}
	else
		qDebug() << "open sucess";

	QFileInfo pic(filename);

	while (!mtlfile.atEnd())
	{
		QByteArray linedate = mtlfile.readLine();

		linedate = linedate.trimmed();

		if (linedate == "")
			continue;

		QList<QByteArray> str_values = linedate.split(' ');
		str_values.removeAll("");
		QString data_type = str_values.takeFirst();

		if (data_type == "newmtl")
		{
			model_3d.mtl_count++;
			Material mtl;
			model_3d.mtls << mtl;
			model_3d.mtls[model_3d.mtl_count - 1].mtlname = str_values[0];
		}

		else if (data_type == "map_Kd")
		{
			model_3d.mtls[model_3d.mtl_count - 1].map_kd_path = pic.path() + "/" + str_values[0];
		}
		else if (data_type == "map_Ks")
		{
			model_3d.mtls[model_3d.mtl_count - 1].map_ks_path = pic.path() + "/" + str_values[0];
		}

		else if (data_type == "Kd")
		{
			model_3d.mtls[model_3d.mtl_count - 1].Kd.setX(str_values[0].toDouble());
			model_3d.mtls[model_3d.mtl_count - 1].Kd.setY(str_values[1].toDouble());
			model_3d.mtls[model_3d.mtl_count - 1].Kd.setZ(str_values[2].toDouble());
		}
		else if (data_type == "Ka")
		{
			model_3d.mtls[model_3d.mtl_count - 1].Ka.setX(str_values[0].toDouble());
			model_3d.mtls[model_3d.mtl_count - 1].Ka.setY(str_values[1].toDouble());
			model_3d.mtls[model_3d.mtl_count - 1].Ka.setZ(str_values[2].toDouble());
		}
		else if (data_type == "Ks")
		{
			model_3d.mtls[model_3d.mtl_count - 1].Ks.setX(str_values[0].toDouble());
			model_3d.mtls[model_3d.mtl_count - 1].Ks.setY(str_values[1].toDouble());
			model_3d.mtls[model_3d.mtl_count - 1].Ks.setZ(str_values[2].toDouble());
		}

		else
			continue;

	}


	mtlfile.close();
		for (int i = 0; i < model_3d.mesh_count; i++)
		{

			for (int j = 0; j < model_3d.mtl_count; j++)
			{

				if (model_3d.meshes[i].mtl.mtlname == model_3d.mtls[j].mtlname)
				{
					model_3d.meshes[i].mtl.Ka = model_3d.mtls[j].Ka;
					model_3d.meshes[i].mtl.Kd = model_3d.mtls[j].Kd;
					model_3d.meshes[i].mtl.Ks = model_3d.mtls[j].Ks;
				}

				if (QFile::exists(model_3d.mtls[j].map_kd_path))
				{
					model_3d.meshes[i].mtl.map_kd_path = model_3d.mtls[j].map_kd_path;
					model_3d.meshes[i].mtl.value_or_diffmap = true;
				}
				else
				{
					model_3d.meshes[i].mtl.value_or_diffmap = false;
				}

				if (QFile::exists(model_3d.mtls[j].map_ks_path))
				{
					model_3d.meshes[i].mtl.map_ks_path = model_3d.mtls[j].map_ks_path;
					model_3d.meshes[i].mtl.value_or_specmap = true;
				}
				else
				{
					model_3d.meshes[i].mtl.value_or_specmap = false;
				}

			}

		}

	return true;
}


void ModelLoder::setFileName(QString filename)
{
	obj_file_path = filename;
}

Model3D ModelLoder::getloadedModel()
{
	return model_3d;
}


void ModelLoder::loadBinocularModel()
{

	//cv::Mat left = cv::imread("left.jpg", cv::IMREAD_GRAYSCALE);
	//cv::Mat right = cv::imread("right.jpg", cv::IMREAD_GRAYSCALE);
	//cv::Mat disp;
	//int mindisparity = 0;
	//int ndisparities = 64;
	//int SADWindowSize = 11;

	////SGBM
	//cv::Ptr<cv::StereoSGBM> sgbm = cv::StereoSGBM::create(mindisparity, ndisparities, SADWindowSize);

	//int P1 = 8 * left.channels() * SADWindowSize* SADWindowSize;
	//int P2 = 32 * left.channels() * SADWindowSize* SADWindowSize;

	//sgbm->setP1(P1);
	//sgbm->setP2(P2);
	//sgbm->setPreFilterCap(5);
	//sgbm->setUniquenessRatio(10);
	//sgbm->setSpeckleRange(2);
	//sgbm->setSpeckleWindowSize(100);
	//sgbm->setDisp12MaxDiff(1);
	////sgbm->setMode(cv::StereoSGBM::MODE_HH);
	//sgbm->compute(left, right, disp);

	//disp.convertTo(disp, CV_32F, 1.0 / 16);                //除以16得到真实视差值
	//cv::Mat disp8U = cv::Mat(disp.rows, disp.cols, CV_8UC1);       //显示
	//normalize(disp, disp8U, 0, 255, cv::NORM_MINMAX, CV_8UC1);


	////cout << disp8U;

	////视差图转为空间点
	//Mesh meshi;
	//model_3d.meshes << meshi;
	//model_3d.mesh_count++;
	//for (int i = 0; i < disp8U.rows; i++)
	//{
	//	for (int j = 0; j < disp8U.cols; j++)
	//	{

	//		uchar disp_value = disp8U.at<uchar>(i, j);


	//		double x = (j / (disp8U.rows / 2.0)) - 1.0;
	//		double y = (i / (disp8U.cols / 2.0)) - 1.0;
	//		//double z = ((255 - disp_value) / 125.0) - 1.0;

	//		double z = 80.0/ disp_value;

	//		model_3d.meshes[0].vps << x << y << z;

	//	}
	//}

	//model_3d.meshes[0].allps << model_3d.meshes[0].vps;
	//
	////qDebug() << model_3d.meshes[0].allps;


	//// 锥体
	//cv::Mat left = cv::imread("steroMatchImage/left10.jpg", cv::IMREAD_GRAYSCALE);
	//cv::Mat right = cv::imread("steroMatchImage/right10.jpg", cv::IMREAD_GRAYSCALE);

	///*立体校正*/
	//cv::Mat cameraMatrix[2], distCoeffs[2];
	//cv::Mat Q;

	//cameraMatrix[0] = cv::Mat((Mat_<double>(3, 3) << 729.0952, 0.0, 0.0,
	//	-0.0867, 729.6895, 0.0,
	//	674.2116, 357.1047, 1.0)).t();

	//distCoeffs[0] = (cv::Mat_<double>(1, 14) << 6.9848610705178477e-02, -3.1045474247642400e-02, -6.5610233314618927e-03, -1.1462353931275203e-04, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.);

	//cameraMatrix[1] = Mat((Mat_<double>(3, 3) << 727.3174, 0.0, 0.0,
	//	-0.2356, 727.4837, 0.0,
	//	634.6699, 342.4619, 1.0)).t();

	//distCoeffs[1] = (cv::Mat_<double>(1, 14) << 9.8797697968316223e-02, -1.2760538137348079e-01,
	//	-8.2947641423894588e-04, 8.3943871476181383e-04, 0., 0., 0., 0.,
	//	0., 0., 0., 0., 0., 0.);

	//Mat R = Mat((Mat_<double>(3, 3) << 1.0, 1.9618e-04, -0.0022,
	//	-1.9726e-04, 1.0, -4.827e-04,
	//	0.0022, 4.83e-04, 1.0)).t();

	//Mat T = (Mat_<double>(3, 1) << -6.1719, -0.0082, 0.0043); //左右相机相互关联的平移向量

	//cv::Mat R1 = (cv::Mat_<double>(3, 3) << 9.9908523344342171e-01, -9.1992112105362279e-03, 4.1762074043410001e-02,
	//	9.3234166309280522e-03, 9.9995267073263494e-01, -2.7803231820107406e-03,
	//	-4.1734520694857037e-02, 3.1671450510257367e-03, 9.9912371555007995e-01);

	//cv::Mat R2 = (cv::Mat_<double>(3, 3) << 9.9913038347395211e-01, -8.7274676713489392e-03, 4.0771413113096469e-02,
	//	8.6212419619687286e-03, 9.9995897066935735e-01, 2.7804972418828274e-03,
	//	-4.0794006989095910e-02, -2.4265790579512678e-03, 9.9916463143360401e-01);

	//cv::Mat Pro1 = (cv::Mat_<double>(3, 4) << 6.1366085042324494e+02, 0., 5.9218479919433594e+02, 0., 0.,
	//	6.1366085042324494e+02, 3.3312004470825195e+02, 0., 0., 0., 1.,
	//	0.);

	//cv::Mat Pro2 = (cv::Mat_<double>(3, 4) << 6.1366085042324494e+02, 0., 5.9218479919433594e+02,
	//	-3.7143416039237313e+03, 0., 6.1366085042324494e+02,
	//	3.3312004470825195e+02, 0., 0., 0., 1., 0.);

	//cv::Rect validRoi[2];
	//cv::Size imageSize = cv::Size(1280, 720);

	////stereoRectify根据内参和畸变系数，计算右相机相对左相机的旋转R和平移矩阵T
	////并将旋转与平移矩阵分解为左右相机各旋转一半的旋转矩阵R1，R2，左右相机在新坐标系下的投影矩阵P1，P2；深度差异映射矩阵Q
	////这里用的是bougust极线校准方法
	//stereoRectify(cameraMatrix[0], distCoeffs[0],
	//	cameraMatrix[1], distCoeffs[1],
	//	imageSize, R, T, R1, R2, Pro1, Pro2, Q,
	//	cv::CALIB_ZERO_DISPARITY, 1, imageSize, &validRoi[0], &validRoi[1]); //validPixROI1  validPixROI2 - 可选的输出参数，Rect型数据。其内部的所有像素都有效

	////计算左右视图的校正查找映射表
	//cv::Mat rmap[2][2];
	//initUndistortRectifyMap(cameraMatrix[0], distCoeffs[0], R1, Pro1, imageSize, CV_16SC2, rmap[0][0], rmap[0][1]);
	//initUndistortRectifyMap(cameraMatrix[1], distCoeffs[1], R2, Pro2, imageSize, CV_16SC2, rmap[1][0], rmap[1][1]);

	////对源图像进行重映射，映射输出为校正后的图像
	//cv::Mat rleftimg, rrightimg;
	//remap(left, rleftimg, rmap[0][0], rmap[0][1], cv::INTER_LINEAR);
	//remap(right, rrightimg, rmap[1][0], rmap[1][1], cv::INTER_LINEAR);

	//cv::Mat sgbmleftimg = rleftimg.clone();
	//cv::Mat sgbmrightimg = rrightimg.clone();

	//cvtColor(rleftimg, rleftimg, cv::COLOR_GRAY2BGR);
	//cvtColor(rrightimg, rrightimg, cv::COLOR_GRAY2BGR);

	////用矩形框起有效区域
	//cv::Rect vroileft(cvRound(validRoi[0].x), cvRound(validRoi[0].y), cvRound(validRoi[0].width), cvRound(validRoi[0].height));
	//rectangle(rleftimg, vroileft, cv::Scalar(0, 0, 255), 3, 8);
	//cv::Rect vroiright(cvRound(validRoi[1].x), cvRound(validRoi[1].y), cvRound(validRoi[1].width), cvRound(validRoi[1].height));
	//rectangle(rrightimg, vroiright, cv::Scalar(0, 0, 255), 3, 8);

	////画出极线进行对比
	//cv::Mat canvas;
	//int w = imageSize.width;
	//int h = imageSize.height;
	//canvas.create(h, w * 2, CV_8UC3);
	//hconcat(rleftimg, rrightimg, canvas);

	//for (int j = 0; j < canvas.rows; j += 16)
	//	line(canvas, cv::Point(0, j), cv::Point(canvas.cols, j), cv::Scalar(0, 255, 0), 1, 8);
	//
	//////显示立体校正结果
	////cv::namedWindow("result_recfy", cv::WINDOW_AUTOSIZE);
	////cv::imshow("result_recfy",canvas);
	////cv::waitKey(0);


	///*立体匹配*/
	//cv::Mat disp, disp8;
	//int mindisparity = 0;
	//int ndisparities = 64;
	//int SADWindowSize = 11;

	////SGBM
	//cv::Ptr<cv::StereoSGBM> sgbm = cv::StereoSGBM::create(mindisparity, ndisparities, SADWindowSize);

	//int P1 = 8 * sgbmleftimg.channels() * SADWindowSize* SADWindowSize;
	//int P2 = 32 * sgbmleftimg.channels() * SADWindowSize* SADWindowSize;

	//sgbm->setP1(P1);
	//sgbm->setP2(P2);
	//sgbm->setMinDisparity(mindisparity);
	//sgbm->setNumDisparities(ndisparities);
	//sgbm->setPreFilterCap(5);
	//sgbm->setUniquenessRatio(10);
	//sgbm->setSpeckleRange(2);
	//sgbm->setSpeckleWindowSize(100);
	//sgbm->setDisp12MaxDiff(1);
	////sgbm->setMode(cv::StereoSGBM::MODE_HH);
	//sgbm->compute(sgbmleftimg, sgbmrightimg, disp);

	//disp.convertTo(disp8, CV_8U, 255 / (ndisparities *16.));

	////显示立体匹配结果
	//cv::Mat disp8U = cv::Mat(disp.rows, disp.cols, CV_8UC1);       //显示
	//cv::normalize(disp, disp8U, 0, 255, cv::NORM_MINMAX, CV_8UC1);
	//cv::namedWindow("result_match", cv::WINDOW_AUTOSIZE);
	//cv::imshow("result_match",disp8U);
	//cv::waitKey(0);


	//cv::Mat xyz;  //三维坐标
	//reprojectImageTo3D(disp, xyz, Q, true); //在实际求距离时，ReprojectTo3D出来的X / W, Y / W, Z / W都要乘以16(也就是W除以16)，才能得到正确的三维坐标信息。
	//xyz = xyz * 16;


	// 重建汽车
	const int imageWidth = 1920;                             //摄像头的分辨率  
	const int imageHeight = 1024;
	Size imageSize = Size(imageWidth, imageHeight);

	Mat rgbImageL, grayImageL;
	Mat rgbImageR, grayImageR;
	Mat rectifyImageL, rectifyImageR;

	Rect validROIL;//图像校正之后，会对图像进行裁剪，这里的validROI就是指裁剪之后的区域  
	Rect validROIR;

	Mat mapLx, mapLy, mapRx, mapRy;     //映射表  
	Mat Rl, Rr, Pl, Pr, Q;              //校正旋转矩阵R，投影矩阵P 重投影矩阵Q
	Mat xyz;              //三维坐标

	Ptr<StereoSGBM> sgbm = StereoSGBM::create(0, 16, 3);

	/*
	事先标定好的相机的参数
	fx 0 cx
	0 fy cy
	0 0  1
	*/

	Mat cameraMatrixL = (Mat_<double>(3, 3) << 4334.09568, 0, 959.50000,//左内参
		0, 4334.09568, 511.50000,
		0, 0, 1.0);
	Mat distCoeffL = (Mat_<double>(5, 1) << 0.0, 0.0, 0.0, 0.0, 0.0);//左畸变系数

	Mat cameraMatrixR = (Mat_<double>(3, 3) << 4489.55770, 0, 801.86552,//右内参
		0, 4507.13412, 530.72579,
		0, 0, 1.0);
	Mat distCoeffR = (Mat_<double>(5, 1) << 0.0, 0.0, 0.0, 0.0, 0.0);//右畸变系数

	Mat T = (Mat_<double>(3, 1) << -518.97666, 01.20629, 9.14632);//左右相机相互关联的平移向量
	Mat rec = (Mat_<double>(3, 1) << 0.04345, -0.05236, -0.01810);//左右相机  相互关联的旋转向量
	Mat R;//左右相机相互关联的旋转矩阵

		/*	立体校正	*/
	Rodrigues(rec, R); //Rodrigues变换      
	stereoRectify(cameraMatrixL, distCoeffL, cameraMatrixR, distCoeffR, imageSize, R, T, Rl, Rr, Pl, Pr, Q, CALIB_ZERO_DISPARITY, 0, imageSize, &validROIL, &validROIR);
	initUndistortRectifyMap(cameraMatrixL, distCoeffL, Rl, Pl, imageSize, CV_16SC2, mapLx, mapLy);
	initUndistortRectifyMap(cameraMatrixR, distCoeffR, Rr, Pr, imageSize, CV_16SC2, mapRx, mapRy);

	/*	读取图片	*/
	rgbImageL = imread("left_cor.bmp", 1);//CV_LOAD_IMAGE_COLOR
	rgbImageR = imread("right_cor.bmp", 1);

	/*	经过remap之后，左右相机的图像已经共面并且行对准了	*/
	remap(rgbImageL, rectifyImageL, mapLx, mapLy, INTER_LINEAR);//INTER_LINEAR
	remap(rgbImageR, rectifyImageR, mapRx, mapRy, INTER_LINEAR);

	///*	把校正结果显示出来*/
	////显示在同一张图上
	//Mat canvas;
	//double sf;
	//int w, h;
	//sf = 700. / MAX(imageSize.width, imageSize.height);
	//w = cvRound(imageSize.width * sf);
	//h = cvRound(imageSize.height * sf);
	//canvas.create(h, w * 2, CV_8UC3);   //注意通道
	//Mat canvasPart = canvas(Rect(w * 0, 0, w, h));                                //得到画布的一部分  
	//resize(rectifyImageL, canvasPart, canvasPart.size(), 0, 0, INTER_AREA);     //把图像缩放到跟canvasPart一样大小  
	//Rect vroiL(cvRound(validROIL.x*sf), cvRound(validROIL.y*sf),                //获得被截取的区域    
	//	cvRound(validROIL.width*sf), cvRound(validROIL.height*sf));
	////rectangle(canvasPart, vroiL, Scalar(0, 0, 255), 3, 8);                      //画上一个矩形  
	//cout << "Painted ImageL" << endl;

	////右图像画到画布上
	//canvasPart = canvas(Rect(w, 0, w, h));                                      //获得画布的另一部分  
	//resize(rectifyImageR, canvasPart, canvasPart.size(), 0, 0, INTER_LINEAR);
	//Rect vroiR(cvRound(validROIR.x * sf), cvRound(validROIR.y*sf),
	//	cvRound(validROIR.width * sf), cvRound(validROIR.height * sf));
	////rectangle(canvasPart, vroiR, Scalar(0, 0, 255), 3, 8);
	//cout << "Painted ImageR" << endl;

	////画上对应的线条
	//for (int i = 0; i < canvas.rows; i += 16)
	//	line(canvas, Point(0, i), Point(canvas.cols, i), Scalar(0, 255, 0), 1, 8);

	//namedWindow("rectified", WINDOW_AUTOSIZE);
	//imshow("rectified", canvas);
	//waitKey(0);


	/*	立体匹配	*/
	sgbm->setPreFilterCap(63);
	int sgbmWinSize = 5;//根据实际情况自己设定
	int NumDisparities = 416;//根据实际情况自己设定
	int UniquenessRatio = 6;//根据实际情况自己设定
	sgbm->setBlockSize(sgbmWinSize);
	int cn = rectifyImageL.channels();

	sgbm->setP1(8 * cn*sgbmWinSize*sgbmWinSize);
	sgbm->setP2(32 * cn*sgbmWinSize*sgbmWinSize);
	sgbm->setMinDisparity(0);
	sgbm->setNumDisparities(NumDisparities);
	sgbm->setUniquenessRatio(UniquenessRatio);
	sgbm->setSpeckleWindowSize(100);
	sgbm->setSpeckleRange(10);
	sgbm->setDisp12MaxDiff(1);
	sgbm->setMode(StereoSGBM::MODE_SGBM);
	Mat disp, disp8, depth;
	sgbm->compute(rectifyImageL, rectifyImageR, disp);

	//深度图
	Mesh meshi;
	model_3d.meshes << meshi;
	model_3d.mesh_count++;
	for (int i = 0; i < disp.rows; i++)
	{
		for (int j = 0; j < disp.cols; j++)
		{
			double x;
			double y;
			double z = disp.at<short int>(i, j);
			if (z > 0)
			{
				z = (1 / z) * 5646 * 0.52;
				//x = (double)i * z / (double)5646 / 1000;
				//y = (double)j * z / (double)5646 / 1000;
				x = (double)i/(double)imageWidth;
				y = (double)j/(double)imageHeight;
			}
			else
			{
				x = 0;
				z = 0;
				y = 0;
			}

			//double x = (double)i/(double)1920;
			//double y = (double)j/(double)1024;


			model_3d.meshes[0].vps << x << y << z;
		}
	}

	////视差图
	//disp.convertTo(disp8, CV_8U, 255 / (NumDisparities *16.));
	//reprojectImageTo3D(disp, xyz, Q, true); //在实际求距离时，ReprojectTo3D出来的X / W, Y / W, Z / W都要乘以16(也就是W除以16)，才能得到正确的三维坐标信息。
	//xyz = xyz * 16;
	////显示视差图
	//namedWindow("disparity", WINDOW_AUTOSIZE);
	//imshow("disparity", disp8);
	//waitKey(0);
	////Mat三维点转opengl
	//Mesh meshi;
	//model_3d.meshes << meshi;
	//model_3d.mesh_count++;
	//for (int i = 0; i < xyz.rows; i++)
	//{
	//	for (int j = 0; j < xyz.cols; j++)
	//	{

	//		double x = xyz.at<cv::Vec3f>(i, j)[0];
	//		double y = xyz.at<cv::Vec3f>(i, j)[1];
	//		double z = xyz.at<cv::Vec3f>(i, j)[2];

	//		model_3d.meshes[0].vps << x << y << z;
	//	}
	//}

	model_3d.meshes[0].allps << model_3d.meshes[0].vps;

}




ModelLoder::~ModelLoder()
{ 
}
