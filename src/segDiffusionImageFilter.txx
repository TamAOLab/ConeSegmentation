/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkGeodesicActiveContourLevelSetImageFilter.txx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __segDiffusionImageFilter_txx
#define __segDiffusionImageFilter_txx

#include "segDiffusionImageFilter.h"

namespace itk {


template< class TInputImage, class TOutputImage>
segDiffusionImageFilter<TInputImage, TOutputImage>
::segDiffusionImageFilter()
{
	m_TimeStep = 1.; // time step 
	m_NumberOfIterations = 10; //number of iterations for level-set propagation
	m_DiffusionWeight = 1.;
}
 
template <class TInputImage, class TOutputImage>
void segDiffusionImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream &os, Indent indent) const
{
	Superclass::PrintSelf(os, indent);
}


template <class TInputImage, class TOutputImage>
void segDiffusionImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
	DiffuseImage();
}

template <class TInputImage, class TOutputImage>
template <class TInput, class TOutput>
typename TOutput::Pointer segDiffusionImageFilter<TInputImage, TOutputImage>
::AllocateImage(typename TInput::Pointer in, typename TOutput::PixelType value)
{
	typename TOutput::Pointer out = TOutput::New();
	out->SetRegions( in->GetLargestPossibleRegion() );
	out->SetSpacing( in->GetSpacing() );
	out->CopyInformation( in );
	out->Allocate();
	out->FillBuffer(value);
	return out;
}

template <class TInputImage, class TOutputImage>
template <class TInput, class TOutput>
typename TOutput::Pointer segDiffusionImageFilter<TInputImage, TOutputImage>
::CopyImage(typename TInput::Pointer input, typename TOutput::Pointer output)
{
	typedef itk::ImageRegionIteratorWithIndex<TInput>	TIterator1;
	typedef itk::ImageRegionIteratorWithIndex<TOutput>	TIterator2;

	TIterator1 it1(input, input->GetLargestPossibleRegion());
	TIterator2 it2(output, output->GetLargestPossibleRegion());

	for (it1.GoToBegin(), it2.GoToBegin(); !it1.IsAtEnd(); ++it1, ++it2)
		it2.Set(it1.Get());

	return output;
}

template <class TInputImage, class TOutputImage>
template <class TInput, class TOutput>
typename TOutput::Pointer segDiffusionImageFilter<TInputImage, TOutputImage>
::CastImage(typename TInput::Pointer in)
{
	typedef itk::CastImageFilter< TInput, TOutput > CastFilterType;
    typename CastFilterType::Pointer castFilter = CastFilterType::New();
	castFilter->SetInput(in);
	castFilter->Update();
	return castFilter->GetOutput();
}

template <class TInputImage, class TOutputImage>
void segDiffusionImageFilter<TInputImage, TOutputImage>
::DiffuseImage()
{
	//first peform mean curvature motion as a starting point
	//step 1: get number of image pixels excluding image boundaries
	InputImagePointer inputPtr = const_cast<InputImageType *>(this->GetInput());
    typename InputImageType::SizeType img_size = inputPtr->GetLargestPossibleRegion().GetSize();
	unsigned int number_of_points = 1;
	for (unsigned int i=0; i<ImageDimension; i++)
		number_of_points *= img_size[i]-2;
		
	//step 2: assign some parameters for AOS computation
	AOSVectorType aos_diagon[ImageDimension];
	AOSVectorType aos_top_diagon[ImageDimension];
	AOSVectorType aos_bottom_diagon[ImageDimension];
	
	AOSVectorType parameter_m;
	AOSVectorType parameter_l;
	AOSVectorType parameter_y;
	
	//we only set two matrix item because AOS computation is perfomed image direction by image direction
	//at each image direction, there are only two neighborhoods
	InternalRealImageIndexType neighbor_index[2];
	InternalRealImageIndexType tmp_neighbor_index[4];
	
	double gradient_values[2], diffusion_values[2];
	double episilon = 0.001;
	double img_dim = ImageDimension;
	unsigned int aos_vector_index, line_index;
	double current_value;
	
	for (unsigned int i=0; i<ImageDimension; i++)
	{
		aos_diagon[i].SetSize(number_of_points);
		aos_top_diagon[i].SetSize(number_of_points-1);
		aos_bottom_diagon[i].SetSize(number_of_points-1);
	}
	parameter_m.SetSize(number_of_points);
	parameter_l.SetSize(number_of_points-1);
	parameter_y.SetSize(number_of_points);
	
	//step 3: AOS computation
	InternalRealImagePointer current_img = AllocateImage<InputImageType, InternalRealImageType>(inputPtr, 0);
	InternalRealImagePointer previous_img = AllocateImage<InputImageType, InternalRealImageType>(inputPtr, 0);
	previous_img = CopyImage<InputImageType, InternalRealImageType>(inputPtr, previous_img);
	
	//AOS computation will only focus on image region without boundary
    InternalRealImageIndexType desiredStart;
	typename InternalRealImageType::SizeType desiredSize;
	for (unsigned int i=0; i<ImageDimension; i++)
	{
		desiredStart[i] = 1;
		desiredSize[i] = inputPtr->GetLargestPossibleRegion().GetSize()[i]-2;
	}
	InternalRealImageRegionType desiredRegion(desiredStart, desiredSize);
	
	for (unsigned int n=0; n<m_NumberOfIterations; n++)
	{
		//(a) initialize aos array
		for (unsigned int j=0; j<ImageDimension; j++)
		{
			for (unsigned int k=0; k<number_of_points-1; k++)
			{
				aos_diagon[j][k] = 0;
				aos_top_diagon[j][k] = 0;
				aos_bottom_diagon[j][k] = 0;
			}
			aos_diagon[j][number_of_points-1] = 0;
		}
		
		//(b) AOS computation	
		//initialize current_levelset_img
		current_img->FillBuffer(0);
		
		InternalRealLinearIteratorType previous_it(previous_img, desiredRegion);
		
		//assign the sparse matrix
		for (unsigned int i=0; i<ImageDimension; i++)
		{
			previous_it.SetDirection(i);
			
			//build aos matrix along different directions
			for (previous_it.GoToBegin(), aos_vector_index = 0; !previous_it.IsAtEnd(); previous_it.NextLine())
			{
				previous_it.GoToBeginOfLine();
				line_index = 0;				
				while ( ! previous_it.IsAtEndOfLine() )
				{
					//at the begining of the current line
					if (line_index == 0)
					{
						for (unsigned int j=0; j<ImageDimension; j++)
							neighbor_index[1][j] = previous_it.GetIndex()[j];
						neighbor_index[1][i] = previous_it.GetIndex()[i] + 1;
						
						gradient_values[1] = 0;
						//compute gradient magnitude along current iteration image direction
						gradient_values[1] += (previous_img->GetPixel(neighbor_index[1]) - previous_it.Get())
							* (previous_img->GetPixel(neighbor_index[1]) - previous_it.Get());
							
						//Accumuate gradient magnitude values along all other image directions
						for (unsigned int m=0; m<ImageDimension; m++)
						{
							if (m == i)
								continue;
								
							for (unsigned int k=0; k<2; k++)
								for (unsigned int j=0; j<ImageDimension; j++)
								{
									tmp_neighbor_index[2*k][j] = previous_it.GetIndex()[j];
									tmp_neighbor_index[2*k+1][j] = neighbor_index[1][j];
								}
							
							tmp_neighbor_index[0][m] = previous_it.GetIndex()[m]-1;
							tmp_neighbor_index[1][m] = previous_it.GetIndex()[m]-1;
							tmp_neighbor_index[2][m] = previous_it.GetIndex()[m]+1;
							tmp_neighbor_index[3][m] = previous_it.GetIndex()[m]+1;
							
							gradient_values[1] += (0.25*(previous_img->GetPixel(tmp_neighbor_index[2])-previous_img->GetPixel(tmp_neighbor_index[0]))
											   +0.25*(previous_img->GetPixel(tmp_neighbor_index[3])-previous_img->GetPixel(tmp_neighbor_index[1])))
											   *(0.25*(previous_img->GetPixel(tmp_neighbor_index[2])-previous_img->GetPixel(tmp_neighbor_index[0]))
											   +0.25*(previous_img->GetPixel(tmp_neighbor_index[3])-previous_img->GetPixel(tmp_neighbor_index[1])));
						}
						
						diffusion_values[1] = 1./pow(gradient_values[1]+episilon*episilon, m_DiffusionWeight/2);
						
						aos_bottom_diagon[i][aos_vector_index] = -img_dim*img_dim*m_TimeStep*diffusion_values[1];
						aos_diagon[i][aos_vector_index] = img_dim+img_dim*img_dim*m_TimeStep*diffusion_values[1];	
					}
					else if (line_index == img_size[i] - 3)
					{
						for (unsigned int j=0; j<ImageDimension; j++)
							neighbor_index[0][j] = previous_it.GetIndex()[j];
						neighbor_index[0][i] = previous_it.GetIndex()[i] - 1;
						
						gradient_values[0] = 0;
						//compute gradient magnitude along current iteration image direction
						gradient_values[0] += (previous_img->GetPixel(neighbor_index[0]) - previous_it.Get())
							* (previous_img->GetPixel(neighbor_index[0]) - previous_it.Get());
							
						//Accumuate gradient magnitude values along all other image directions
						for (unsigned int m=0; m<ImageDimension; m++)
						{
							if (m == i)
								continue;
								
							for (unsigned int k=0; k<2; k++)
								for (unsigned int j=0; j<ImageDimension; j++)
								{
									tmp_neighbor_index[2*k][j] = previous_it.GetIndex()[j];
									tmp_neighbor_index[2*k+1][j] = neighbor_index[0][j];
								}
							
							tmp_neighbor_index[0][m] = previous_it.GetIndex()[m]-1;
							tmp_neighbor_index[1][m] = previous_it.GetIndex()[m]-1;
							tmp_neighbor_index[2][m] = previous_it.GetIndex()[m]+1;
							tmp_neighbor_index[3][m] = previous_it.GetIndex()[m]+1;
							
							gradient_values[0] += (0.25*(previous_img->GetPixel(tmp_neighbor_index[2])-previous_img->GetPixel(tmp_neighbor_index[0]))
											   +0.25*(previous_img->GetPixel(tmp_neighbor_index[3])-previous_img->GetPixel(tmp_neighbor_index[1])))
											   *(0.25*(previous_img->GetPixel(tmp_neighbor_index[2])-previous_img->GetPixel(tmp_neighbor_index[0]))
											   +0.25*(previous_img->GetPixel(tmp_neighbor_index[3])-previous_img->GetPixel(tmp_neighbor_index[1])));
						}
						
						diffusion_values[0] = 1./pow(gradient_values[0]+episilon*episilon, m_DiffusionWeight/2);
						
						aos_top_diagon[i][aos_vector_index-1] = -img_dim*img_dim*m_TimeStep*diffusion_values[0];
						aos_diagon[i][aos_vector_index] = img_dim+img_dim*img_dim*m_TimeStep*diffusion_values[0];	
						/*if (isNaN<double>(aos_diagon[i][aos_vector_index]))
							cout << "found here" << endl;		*/
					}
					else
					{
						for (unsigned int j=0; j<ImageDimension; j++)
						{
							neighbor_index[0][j] = previous_it.GetIndex()[j];
							neighbor_index[1][j] = previous_it.GetIndex()[j];
						}
						neighbor_index[0][i] = previous_it.GetIndex()[i] - 1;
						neighbor_index[1][i] = previous_it.GetIndex()[i] + 1;
						
						gradient_values[0] = 0;
						gradient_values[1] = 0;
						//compute gradient magnitude along current iteration image direction
						gradient_values[0] += (previous_img->GetPixel(neighbor_index[0]) - previous_it.Get())
							* (previous_img->GetPixel(neighbor_index[0]) - previous_it.Get());
						gradient_values[1] += (previous_img->GetPixel(neighbor_index[1]) - previous_it.Get())
							* (previous_img->GetPixel(neighbor_index[1]) - previous_it.Get());
						
						//Accumuate gradient magnitude values along all other image directions
						for (unsigned int m=0; m<ImageDimension; m++)
						{
							if (m == i)
								continue;
								
							for (unsigned int k=0; k<2; k++)
								for (unsigned int j=0; j<ImageDimension; j++)
								{
									tmp_neighbor_index[2*k][j] = previous_it.GetIndex()[j];
									tmp_neighbor_index[2*k+1][j] = neighbor_index[0][j];
								}
							
							tmp_neighbor_index[0][m] = previous_it.GetIndex()[m]-1;
							tmp_neighbor_index[1][m] = previous_it.GetIndex()[m]-1;
							tmp_neighbor_index[2][m] = previous_it.GetIndex()[m]+1;
							tmp_neighbor_index[3][m] = previous_it.GetIndex()[m]+1;
							
							gradient_values[0] += (0.25*(previous_img->GetPixel(tmp_neighbor_index[2])-previous_img->GetPixel(tmp_neighbor_index[0]))
											   +0.25*(previous_img->GetPixel(tmp_neighbor_index[3])-previous_img->GetPixel(tmp_neighbor_index[1])))
											   *(0.25*(previous_img->GetPixel(tmp_neighbor_index[2])-previous_img->GetPixel(tmp_neighbor_index[0]))
											   +0.25*(previous_img->GetPixel(tmp_neighbor_index[3])-previous_img->GetPixel(tmp_neighbor_index[1])));
						}
						diffusion_values[0] = 1./pow(gradient_values[0]+episilon*episilon, m_DiffusionWeight/2);
						
						for (unsigned int m=0; m<ImageDimension; m++)
						{
							if (m == i)
								continue;
								
							for (unsigned int k=0; k<2; k++)
								for (unsigned int j=0; j<ImageDimension; j++)
								{
									tmp_neighbor_index[2*k][j] = previous_it.GetIndex()[j];
									tmp_neighbor_index[2*k+1][j] = neighbor_index[1][j];
								}
							
							tmp_neighbor_index[0][m] = previous_it.GetIndex()[m]-1;
							tmp_neighbor_index[1][m] = previous_it.GetIndex()[m]-1;
							tmp_neighbor_index[2][m] = previous_it.GetIndex()[m]+1;
							tmp_neighbor_index[3][m] = previous_it.GetIndex()[m]+1;
							
							gradient_values[1] += (0.25*(previous_img->GetPixel(tmp_neighbor_index[2])-previous_img->GetPixel(tmp_neighbor_index[0]))
											   +0.25*(previous_img->GetPixel(tmp_neighbor_index[3])-previous_img->GetPixel(tmp_neighbor_index[1])))
											   *(0.25*(previous_img->GetPixel(tmp_neighbor_index[2])-previous_img->GetPixel(tmp_neighbor_index[0]))
											   +0.25*(previous_img->GetPixel(tmp_neighbor_index[3])-previous_img->GetPixel(tmp_neighbor_index[1])));
						}
						diffusion_values[1] = 1./pow(gradient_values[1]+episilon*episilon, m_DiffusionWeight/2);
						
						aos_top_diagon[i][aos_vector_index-1] = -img_dim*img_dim*m_TimeStep*diffusion_values[0];
						aos_bottom_diagon[i][aos_vector_index] = -img_dim*img_dim*m_TimeStep*diffusion_values[1];
						aos_diagon[i][aos_vector_index] = img_dim+img_dim*img_dim*m_TimeStep*(diffusion_values[0]+diffusion_values[1]);
					}
					++previous_it;
					line_index++;
					aos_vector_index++;
				}
			}
		}
			
		//(C) solve the tridiagonal linear system using Thomas algorithm
		for (unsigned int i=0; i<ImageDimension; i++)
		{
			//LR decomposition
			parameter_m[0] = aos_diagon[i][0];
			for (unsigned int k=0; k<number_of_points-1; k++)
			{
				parameter_l[k] = aos_bottom_diagon[i][k]/parameter_m[k];
				parameter_m[k+1] = aos_diagon[i][k+1]-parameter_l[k]*aos_top_diagon[i][k];
			}
			
			previous_it.SetDirection(i);
			
			//forward substitution
			for (previous_it.GoToBegin(), aos_vector_index = 0; !previous_it.IsAtEnd(); previous_it.NextLine())
			{		
				previous_it.GoToBeginOfLine();
				while ( ! previous_it.IsAtEndOfLine() )
				{	
					if (aos_vector_index == 0)
					{
						parameter_y[aos_vector_index] = previous_it.Get();
					}
					else
					{
						parameter_y[aos_vector_index] = previous_it.Get() - parameter_l[aos_vector_index-1]*parameter_y[aos_vector_index-1];
					}
						
					++previous_it;
					aos_vector_index ++;
				}
			}
			
			//backward substitution
			for (previous_it.GoToReverseBegin(), aos_vector_index = number_of_points-1; !previous_it.IsAtReverseEnd(); previous_it.PreviousLine())
			{		
				previous_it.GoToReverseBeginOfLine();
				while ( ! previous_it.IsAtReverseEndOfLine() )
				{
					//cout << parameter_y[aos_vector_index] << ", " << parameter_m[aos_vector_index] << endl;
					if (aos_vector_index == number_of_points-1)
					{
						current_value = parameter_y[aos_vector_index]/parameter_m[aos_vector_index];
					}
					else
					{
						current_value = (parameter_y[aos_vector_index]-aos_top_diagon[i][aos_vector_index]*current_value)/parameter_m[aos_vector_index];
					}
					current_img->SetPixel(previous_it.GetIndex(), current_img->GetPixel(previous_it.GetIndex())+current_value);
					
					--previous_it;
					aos_vector_index --;
				}
			}
		}
		
		for (unsigned int i=0; i<ImageDimension; i++)
		{
			//One image side
            InternalRealImageIndexType boundaryStart, boundaryStart1;
            typename InternalRealImageType::SizeType boundarySize, boundarySize1;
			for (unsigned int j=0; j<ImageDimension; j++)
			{
				boundaryStart[j] = 0;
				boundaryStart1[j] = 0;
				boundarySize[j] = inputPtr->GetLargestPossibleRegion().GetSize()[j];
				boundarySize1[j] = inputPtr->GetLargestPossibleRegion().GetSize()[j];
			
				if (i == j)
				{
					boundaryStart1[j] = 1;
					boundarySize[j] = boundarySize1[j] = 1;
				}
			}
			
			InternalRealImageRegionType boundaryRegion(boundaryStart, boundarySize);		
			InternalRealImageRegionType boundaryRegion1(boundaryStart1, boundarySize1);
			InternalRealImageIteratorType boundary_it(current_img, boundaryRegion);
			InternalRealImageIteratorType boundary_it1(current_img, boundaryRegion1);
				
			for (boundary_it.GoToBegin(), boundary_it1.GoToBegin(); !boundary_it.IsAtEnd(); ++boundary_it, ++boundary_it1)
				boundary_it.Set(boundary_it1.Get());
			
			//the other image side
			for (unsigned int j=0; j<ImageDimension; j++)
			{
				boundaryStart[j] = 0;
				boundaryStart1[j] = 0;
				boundarySize[j] = inputPtr->GetLargestPossibleRegion().GetSize()[j];
				boundarySize1[j] = inputPtr->GetLargestPossibleRegion().GetSize()[j];
				if (i == j)
				{
					boundaryStart[j] = inputPtr->GetLargestPossibleRegion().GetSize()[j]-2;
					boundaryStart1[j] = inputPtr->GetLargestPossibleRegion().GetSize()[j]-3;
					boundarySize[j] = boundarySize1[j] = 1;
				}
			}
			boundaryRegion.SetIndex(boundaryStart);
			boundaryRegion.SetSize(boundarySize);	
			boundaryRegion1.SetIndex(boundaryStart1);	
			boundaryRegion1.SetSize(boundarySize1);
			boundary_it = InternalRealImageIteratorType(current_img, boundaryRegion);
			boundary_it1 = InternalRealImageIteratorType(current_img, boundaryRegion1);
				
			for (boundary_it.GoToBegin(), boundary_it1.GoToBegin(); !boundary_it.IsAtEnd(); ++boundary_it, ++boundary_it1)
				boundary_it.Set(boundary_it1.Get());
		}
			
		//(d): copy current image to the m_levelset
		previous_img = CopyImage<InternalRealImageType, InternalRealImageType>(current_img, previous_img);
		//WriteImage<InternalRealImageType>("current_img.hdr", current_img);
	}
	
	OutputImagePointer outputPtr = CastImage<InternalRealImageType, TOutputImage>(current_img);
	//WriteImage<TOutputImage>("output1.hdr", outputPtr);
	this->GraftOutput(outputPtr);
	//WriteImage<TOutputImage>("output2.hdr", outputPtr);
}

}// end namespace itk

#endif
