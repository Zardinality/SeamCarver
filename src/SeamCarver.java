import java.awt.Color;

import edu.princeton.cs.algs4.Picture;
import edu.princeton.cs.algs4.Stack;

public class SeamCarver {
	private Picture target_pic;
	private int height;
	private int width;
	private int[][] red_pixel_arr;
	private int[][] green_pixel_arr;
	private int[][] blue_pixel_arr;
	private int[][] edgeTo;
	private double[][] distance;
	private double[][] energy;
	private boolean isvert;

	public SeamCarver(Picture picture) {
		target_pic = new Picture(picture);
		height = target_pic.height();
		width = target_pic.width();
		isvert = true; // is the picture vertical
		red_pixel_arr = new int[width][height];
		green_pixel_arr = new int[width][height];
		blue_pixel_arr = new int[width][height];
		energy = new double[width][height];
		pic2arr(target_pic);
		energyInit();
	} // create a seam carver object based on the given picture

	public Picture picture() {
		if (!isvert)
			changeState();
		target_pic = new Picture(width, height);
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				Color temp = new Color(red_pixel_arr[i][j],
						green_pixel_arr[i][j], blue_pixel_arr[i][j]);
				target_pic.set(i, j, temp);
			}
		}
		return target_pic;
	} // current picture

	public int width() {
		return width;
	}

	public int height() {
		return height;
	}

	private int widthInside() {
		if (isvert)
			return width;
		return height;
	} // width of current picture

	private int heightInside() {
		if (isvert)
			return height;
		return width;
	} // height of current picture

	public double energy(int x, int y) {
		// if(!isvert){
		// int temp = y;
		// y = x;
		// x = temp;
		// }

		if (x < 0 || y < 0 || x >= width || y >= height)
			throw new IllegalArgumentException("Bad pixel coordinates: x = "
					+ x + "; y = " + y);
		if (x == 0 || y == 0 || x == width - 1 || y == height - 1)
			return 1000;
		return Math.sqrt(deltaColor(x - 1, y, x + 1, y)
				+ deltaColor(x, y - 1, x, y + 1));
	} // energy of pixel at column x and row y

	private double energyinside(int x, int y) {

		if (x < 0 || y < 0 || x >= widthInside() || y >= heightInside())
			throw new IllegalArgumentException("Bad pixel coordinates: x = "
					+ x + "; y = " + y);
		if (x == 0 || y == 0 || x == widthInside() - 1
				|| y == heightInside() - 1)
			return 1000;
		return Math.sqrt(deltaColorInside(x - 1, y, x + 1, y)
				+ deltaColorInside(x, y - 1, x, y + 1));
	} // energy of pixel at column x and row y

	public int[] findHorizontalSeam() {
		if (!isvert)
			changeState();
		int[] result = new int[widthInside()];
		edgeTo = new int[widthInside()][heightInside()];
		distance = new double[widthInside()][heightInside()];
		for (int i = 0; i < heightInside(); i++) {
			distance[0][i] = energy[0][i];
			edgeTo[0][i] = i;
		}
		for (int i = 1; i < widthInside(); i++) {
			for (int j = 0; j < heightInside(); j++)
				relax(i, j);
		}

		Stack<Integer> temp = auxiliary();
		for (int i = 0; i < widthInside(); i++)
			result[i] = temp.pop();
		return result;
	} // sequence of indices for horizontal seam

	public int[] findVerticalSeam() {
		if (isvert)
			changeState();
		int[] result = new int[widthInside()];
		edgeTo = new int[widthInside()][heightInside()];
		distance = new double[widthInside()][heightInside()];
		for (int i = 0; i < heightInside(); i++) {
			distance[0][i] = energy[0][i];
			edgeTo[0][i] = i;
		}
		for (int i = 1; i < widthInside(); i++) {
			for (int j = 0; j < heightInside(); j++)
				relax(i, j);
		}
		Stack<Integer> temp = auxiliary();
		for (int i = 0; i < widthInside(); i++)
			result[i] = temp.pop();
		return result;
	} // sequence of indices for vertical seam

	public void removeHorizontalSeam(int[] seam) {
		
		if (!isvert)
			changeState();

		int[][] new_red_pixel_arr = new int[widthInside()][heightInside()];
		int[][] new_green_pixel_arr = new int[widthInside()][heightInside()];
		int[][] new_blue_pixel_arr = new int[widthInside()][heightInside()];
		double[][] new_energy = new double[widthInside()][heightInside()];
		height--;
		for (int i = 0; i < widthInside(); i++) {
			new_red_pixel_arr[i] = deleteATerm(red_pixel_arr[i], seam[i]);
			new_green_pixel_arr[i] = deleteATerm(green_pixel_arr[i], seam[i]);
			new_blue_pixel_arr[i] = deleteATerm(blue_pixel_arr[i], seam[i]);
			new_energy[i] = deleteATerm(energy[i], seam[i]);

		}
		red_pixel_arr = new_red_pixel_arr;
		blue_pixel_arr = new_blue_pixel_arr;
		green_pixel_arr = new_green_pixel_arr;

		for (int i = 1; i < widthInside(); i++) {

			if(seam[i] > 0 && seam[i] < heightInside()){
			   new_energy[i][seam[i]] = energyinside(i, seam[i]);
			   new_energy[i][seam[i] - 1] = energyinside(i, seam[i] - 1);
			}
			else if(seam[i] != 0){
				throw new IllegalArgumentException("Bad seam" + seam[i]);
			}
			new_energy[i][seam[i]] = energyinside(i, seam[i]);
		}
		energy = new_energy;
		if (energy.length != widthInside()
				|| energy[0].length != heightInside())
			throw new IllegalArgumentException("Bad pixel coordinates: dim1 = "
					+ widthInside() + "; dim2 = " + heightInside()
					+ " energydim" + energy.length + " " + energy[0].length);
	} // remove horizontal seam from current picture

	public void removeVerticalSeam(int[] seam) {
		
		if (isvert)
			changeState();

		int[][] new_red_pixel_arr = new int[widthInside()][heightInside()];
		int[][] new_green_pixel_arr = new int[widthInside()][heightInside()];
		int[][] new_blue_pixel_arr = new int[widthInside()][heightInside()];
		double[][] new_energy = new double[widthInside()][heightInside()];
		width--;
		for (int i = 0; i < widthInside(); i++) {
			new_red_pixel_arr[i] = deleteATerm(red_pixel_arr[i], seam[i]);
			new_green_pixel_arr[i] = deleteATerm(green_pixel_arr[i], seam[i]);
			new_blue_pixel_arr[i] = deleteATerm(blue_pixel_arr[i], seam[i]);
			new_energy[i] = deleteATerm(energy[i], seam[i]);

		}
		red_pixel_arr = new_red_pixel_arr;
		blue_pixel_arr = new_blue_pixel_arr;
		green_pixel_arr = new_green_pixel_arr;
		for (int i = 1; i < widthInside(); i++) {
			// if(seam[i] >= heightInside()-2){
			// new_energy[i][seam[i]] = 1000;
			// new_energy[i][seam[i] - 1] = energy(i, seam[i] - 1);
			// continue;
			// }
			// if(seam[i] <= 2){
			// new_energy[i][seam[i]] = energy(i, seam[i]);
			// new_energy[i][seam[i] - 1] = 1000;
			// continue;
			// }
			if(seam[i] > 0 && seam[i] < heightInside()){
			   new_energy[i][seam[i]] = energyinside(i, seam[i]);
			   new_energy[i][seam[i] - 1] = energyinside(i, seam[i] - 1);
			}
			else if(seam[i] != 0){
				throw new IllegalArgumentException("Bad seam" + seam[i]);
			}
			new_energy[i][seam[i]] = energyinside(i, seam[i]);
		}
		energy = new_energy;
		if (energy.length != widthInside()
				|| energy[0].length != heightInside())
			throw new IllegalArgumentException("Bad pixel coordinates: dim1 = "
					+ widthInside() + "; dim2 = " + heightInside()
					+ " energydim" + energy.length + " " + energy[0].length);
	} // remove vertical seam from current picture

	private void pic2arr(Picture picture) {
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				Color color = target_pic.get(i, j);
				red_pixel_arr[i][j] = color.getRed();
				green_pixel_arr[i][j] = color.getGreen();
				blue_pixel_arr[i][j] = color.getBlue();
			}
		}
	}

	private void energyInit() {
		for (int i = 0; i < widthInside(); i++) {
			energy[i][0] = 1000;
			energy[i][heightInside() - 1] = 1000;
		}
		for (int j = 0; j < heightInside(); j++) {
			energy[0][j] = 1000;
			energy[widthInside() - 1][j] = 1000;
		}
		for (int i = 1; i < widthInside() - 1; i++) {
			for (int j = 1; j < heightInside() - 1; j++)
				energy[i][j] = energyinside(i, j);
		}
	}

	private double deltaColor(int i1, int j1, int i2, int j2) {
		if (!isvert) {
			int temp = i1;
			i1 = j1;
			j1 = temp;
		}
		if (!isvert) {
			int temp = i2;
			i2 = j2;
			j2 = temp;
		}
		int fr = red_pixel_arr[i1][j1];
		int fg = green_pixel_arr[i1][j1];
		int fb = blue_pixel_arr[i1][j1];
		int sr = red_pixel_arr[i2][j2];
		int sg = green_pixel_arr[i2][j2];
		int sb = blue_pixel_arr[i2][j2];
		double jerk = Math.pow(fr - sr, 2) + Math.pow(fg - sg, 2)
				+ Math.pow(fb - sb, 2);
		return jerk;
	}

	private double deltaColorInside(int i1, int j1, int i2, int j2) {
		// if(!isvert){
		// int temp = i1;
		// i1 = j1;
		// j1 = temp;
		// }
		// if(!isvert){
		// int temp = i2;
		// i2 = j2;
		// j2 = temp;
		// }
		int fr = red_pixel_arr[i1][j1];
		int fg = green_pixel_arr[i1][j1];
		int fb = blue_pixel_arr[i1][j1];
		int sr = red_pixel_arr[i2][j2];
		int sg = green_pixel_arr[i2][j2];
		int sb = blue_pixel_arr[i2][j2];
		double jerk = (fr - sr) * (fr - sr) + (fg - sg) * (fg - sg)
				+ (fb - sb) * (fb - sb);
		return jerk;
	}

	private void relax(int i, int j) {
		int min_index = 0;
		double min_dis = Integer.MAX_VALUE;
		if (i < 0 || j < 0 || i > widthInside() - 1 || j > heightInside() - 1) {
			throw new IllegalArgumentException("Bad pixel coordinates: i = "
					+ i + "; j = " + j);
		}
		if (i == 0 || j == 0 || i == widthInside() - 1
				|| j == heightInside() - 1) {
			edgeTo[i][j] = j;
			if (i == 0)
				distance[i][j] = 1000;
			else
				distance[i][j] = distance[i-1][j] + 1000;
			return;
		}
		for (int t = j - 1; t < j + 2; t++) {
			if (distance[i-1][t] < min_dis) {
				min_dis = distance[i-1][t];
				min_index = t;
			}
		}
		edgeTo[i][j] = min_index;
		if (energy.length != widthInside()
				|| energy[0].length != heightInside())
			throw new IllegalArgumentException("Bad pixel coordinates: dim1 = "
					+ widthInside() + "; dim2 = " + heightInside()
					+ " energydim" + energy.length + " " + energy[0].length);
		distance[i][j] = min_dis + energy[i][j];
	}

	private void changeState() {
		red_pixel_arr = transpose(red_pixel_arr);
		green_pixel_arr = transpose(green_pixel_arr);
		blue_pixel_arr = transpose(blue_pixel_arr);
		energy = transpose(energy);
		isvert = !isvert;
		// energy = new double[widthInside()][heightInside()];
		// energyInit();
	}

	private double[][] transpose(double[][] arr) {
		// TODO Auto-generated method stub
		int dim1 = arr.length;
		int dim2 = arr[0].length;
		double[][] result = new double[dim2][dim1];
		for (int i = 0; i < dim1; i++)
			for (int j = 0; j < dim2; j++)
				result[j][i] = arr[i][j];
		return result;
	}

	private int[][] transpose(int[][] arr) {
		int dim1 = arr.length;
		int dim2 = arr[0].length;
		int[][] result = new int[dim2][dim1];
		for (int i = 0; i < dim1; i++)
			for (int j = 0; j < dim2; j++)
				result[j][i] = arr[i][j];
		return result;
	}

	private Stack<Integer> auxiliary() {
		double min_dis = Integer.MAX_VALUE;
		int min_index = 0;
		for (int i = 0; i < heightInside(); i++) {
			if (distance[widthInside() - 1][i] < min_dis) {
				min_dis = distance[widthInside() - 1][i];
				min_index = i;
			}
		}
		Stack<Integer> result = new Stack<Integer>();
		result.push(min_index);
		result.push(min_index);
		for (int j = widthInside() - 1; j > 0; j--) {
			min_index = edgeTo[j][min_index];
			result.push(min_index);
		}	
		return result;
	}

	private int[] deleteATerm(int[] temp, int delete_pos) {
		int[] result = new int[temp.length - 1];
		System.arraycopy(temp, 0, result, 0, delete_pos);
		System.arraycopy(temp, delete_pos + 1, result, delete_pos, temp.length
				- delete_pos - 1);
		return result;
	}

	private double[] deleteATerm(double[] temp, int delete_pos) {
		double[] result = new double[temp.length - 1];
		System.arraycopy(temp, 0, result, 0, delete_pos);
		System.arraycopy(temp, delete_pos + 1, result, delete_pos, temp.length
				- delete_pos - 1);
		return result;
	}
}