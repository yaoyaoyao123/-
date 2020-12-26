#include "ModelLoder.h"
#include<QFileInfo>

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


	//Mat img_disp, disp8;
	//Mat img_left = imread("./left01.jpg", 0);
	//Mat img_right = imread("./right01.jpg", 0);

	//resize(img_left, img_left, Size(), 0.2, 0.2);
	//resize(img_right, img_right, Size(), 0.2, 0.2);

	//Size imgSize = img_left.size();
	//int numberOfDisparities = ((imgSize.width / 8) + 15) & -16;
	//cv::Ptr<cv::StereoSGBM> sgbm = cv::StereoSGBM::create(0, 16, 3);
	//sgbm->setPreFilterCap(32);
	//int SADWindowSize = 9;
	//int sgbmWinSize = SADWindowSize > 0 ? SADWindowSize : 3;
	//sgbm->setBlockSize(sgbmWinSize);
	//int cn = img_left.channels();
	//sgbm->setP1(8 * cn*sgbmWinSize*sgbmWinSize);
	//sgbm->setP2(32 * cn*sgbmWinSize*sgbmWinSize);
	//sgbm->setMinDisparity(0);
	//sgbm->setNumDisparities(numberOfDisparities);
	//sgbm->setUniquenessRatio(10);
	//sgbm->setSpeckleWindowSize(100);
	//sgbm->setSpeckleRange(32);
	//sgbm->setDisp12MaxDiff(1);
	//sgbm->setMode(cv::StereoSGBM::MODE_SGBM);
	//sgbm->compute(img_left, img_right, img_disp);
	//// 显示视差图
	//img_disp.convertTo(disp8, CV_8U, 255 / (numberOfDisparities *16.));
	//imshow("disp", disp8);
	//waitKey(0);

	//Mat depth = imread("disp0.png", IMREAD_UNCHANGED);// depth信息
	//// 相机内参
	//const double u0 = 1288.147;
	//const double v0 = 973.571;
	//const double fx = 4152.073;
	//const double fy = 4152.073;
	//const double baseline = 1.76252; // 
	//const double doffs = 213.084;// 两个相机主点在x方向上的差距, doffs = cx1 - cx0

	////点云
	//Mesh meshi;
	//model_3d.meshes << meshi;
	//model_3d.mesh_count++;
	//for (int row = 0; row < depth.rows; row++)
	//{
	//	for (int col = 0; col < depth.cols; col++)
	//	{
	//		ushort d = depth.ptr<ushort>(row)[col];

	//		//if (d == 0)
	//		//	continue;
	//		// depth			
	//		double z = (fx * baseline / (d + doffs)); // z = baseline * f / (d + doffs)
	//		double x = (col - u0)*z / fx;// Xw向右，Yw向下为正
	//		double y = (row - v0)*z / fy;

	//		y = -y; // 为便于显示，绕x轴三维旋转180°
	//		z = -z;
	//		model_3d.meshes[0].vps << x << y << z;
	//		//model_3d.meshes[0].tps << col / 2872.0 << row / 1954.0;
	//	}
	//}

	////贴纹理
	//model_3d.mtl_count++;
	//model_3d.meshes[0].mtl.map_kd_path = "./dispc1.png";
	//model_3d.meshes[0].mtl.value_or_diffmap = true;
	//model_3d.meshes[0].mtl.value_or_specmap = false;
	

	////点云
	//Mesh meshi;
	//model_3d.meshes << meshi;
	//model_3d.mesh_count++;
	//for (int i = 0; i < img_disp.rows; i++)
	//{
	//	for (int j = 0; j < img_disp.cols; j++)
	//	{
	//		//double x;
	//		//double y;
	//		//double z;
	//		//double disp = img_disp.at<short int>(i, j);
	//		//if (disp > 0)
	//		//{
	//		//	z = (1 / disp) * 720 * 6.0;
	//		//	x = (double)i * z / (double)5646;
	//		//	y = (double)j * z / (double)5646;
	//		//	//x = (double)i/(double)imgSize.width;
	//		//	//y = (double)j/(double)imgSize.height;
	//		//}
	//		//else
	//		//{
	//		//	x = (double)i / (double)imgSize.width;
	//		//	y = (double)j / (double)imgSize.height;
	//		//	z = 0;
	//		//}

	//		double x;
	//		double y;
	//		double z;
	//		// 获取深度图中(m,n)处的值
	//		double d = img_disp.ptr<ushort>(i)[j];
	//		z = (720.0 * 60.0) / d;
	//		x = (i* z) / 720.0;
	//		y = (j* z) / 720.0;

	//		model_3d.meshes[0].vps << x << y << z;

	//	}
	//}

	model_3d.meshes[0].allps << model_3d.meshes[0].tps << model_3d.meshes[0].vps;

}




ModelLoder::~ModelLoder()
{ 
}
