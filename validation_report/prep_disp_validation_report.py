#!/usr/bin/env python3

"""Automated routine to prepare a DISP Cal/Val report"""

import argparse
import os
import shutil
from pathlib import Path

import fitz  # PyMuPDF
import matplotlib.pyplot as plt
import pandas as pd
import PyPDF2
from PIL import Image as PIL_Image
from tqdm import tqdm
from matplotlib.colors import to_rgb  # Convert color names & hex to RGB
from wand.image import Image

def create_parser():
    """ Initiate input parser"""

    parser = argparse.ArgumentParser(description='Create DISP '
                                     'validation report.')
    parser.add_argument('-p', '--pdir', dest='pdir', type=str,
                        help='Specify directory containing pngs')
    parser.add_argument('-o', '--outdir', dest='outdir', type=str,
                        default='./', help='Specify output directory')
    parser.add_argument('-f', '--fname', dest='fname', type=str,
                        default='output.pdf', help='Specify output filename')

    return parser


def cmd_line_parse(iargs=None):
    """ Wrap around input parser"""
    parser = create_parser()
    inputs = parser.parse_args(args=iargs)

    return inputs


def main():
    calval_inps = cmd_line_parse()
    prep_disp_report(calval_inps)


def merge_pdf_pages(input_pdf,
                    output_pdf,
                    width,
                    height,
                    resolution,
                    top=0,
                    left=0):
    """ Merge two PDF pages into one"""
    # Create a new image (blank page)
    with Image(width=width,
               height=(height+top)*2,
               resolution=resolution) as new_page:
        # Read the first page of the input PDF
        with Image(filename=f'{input_pdf}[0]') as page1:
            # Set the dimensions for the new page based on the first page
            new_page.width = page1.width
            new_page.height = page1.height * 2  # Concat 2 pages vertically

            # Composite the first page onto the new page
            new_page.composite(page1, top=top, left=left)

            # Read the second page of the input PDF
            with Image(filename=f'{input_pdf}[1]') as page2:
                # Composite the second page below the first page
                new_page.composite(page2, top=page1.height+top, left=0)

        # add border (color, width, and height)
        new_page.border('White', 12, 0)
        # Save the concatenated PDF
        new_page.save(filename=output_pdf)

    return


def concatenate_pdf_pages(input_pdf,
                          output_pdf):
    """ Consolidate two PDFs into one"""
    # Make copy of input
    input_pdf1 = str(output_pdf).replace('.pdf', '_temp.pdf')
    Path(output_pdf).replace(input_pdf1)

    # Open the input PDF files
    with open(input_pdf1, 'rb') as file1, open(input_pdf, 'rb') as file2:
        # Create PDF reader objects
        pdf_reader1 = PyPDF2.PdfReader(file1)
        pdf_reader2 = PyPDF2.PdfReader(file2)

        # Create a PDF writer object for the output file
        pdf_writer = PyPDF2.PdfWriter()

        # Add all pages from the first PDF
        for page_num in range(len(pdf_reader1.pages)):
            page = pdf_reader1.pages[page_num]
            pdf_writer.add_page(page)

        # Add all pages from the second PDF
        for page_num in range(len(pdf_reader2.pages)):
            page = pdf_reader2.pages[page_num]
            pdf_writer.add_page(page)

        # Save the concatenated PDF to the output file
        with open(output_pdf, 'wb') as output_file:
            pdf_writer.write(output_file)

    # delete temp file
    Path(input_pdf1).unlink()

    return


def add_text(output_pdf, header_txt, font_size=24, font_color="black",
             bck_color="#CCCCCC", pos_x=40, pos_y=60, padding=5):
    """
    Add title text to a PDF with a background fill color.

    Args:
        output_pdf (str): Path to the PDF file.
        header_txt (str): Text to add.
        font_size (int): Font size.
        font_color (str): Named or hex color for text.
        bck_color (str): Named or hex color for background.
        pos_x (float): X-coordinate of text.
        pos_y (float): Y-coordinate of text.
        padding (float): Extra space around text for background.
    """
    def convert_color(color):
        """Convert named/hex colors to PyMuPDF's (R, G, B) format."""
        return to_rgb(color)  # Converts "black" or "#CCCCCC" to (R, G, B)

    temp_pdf = str(output_pdf)[:-4] + '_temp.pdf'

    doc = fitz.open(output_pdf)  # Open the PDF
    page = doc[0]  # Select the first page

    # Convert colors
    font_rgb = convert_color(font_color)
    bck_rgb = convert_color(bck_color)

    # Measure text size
    text_width = fitz.get_text_length(header_txt, fontsize=font_size)
    text_height = font_size* 1.2  # Approximate text height

    # Define background rectangle
    # Adjust rectangle position (shift upward)
    rect = fitz.Rect(
        pos_x - padding,
        pos_y - text_height + padding,  # Adjust upward
        pos_x + text_width + padding,
        pos_y + padding
    )

    # Draw filled background
    page.draw_rect(rect, color=bck_rgb, fill=bck_rgb)

    # Add text on top of the background
    page.insert_text((pos_x, pos_y), header_txt, fontsize=font_size,
        color=font_rgb)

    # Save the modified PDF
    doc.save(temp_pdf)
    doc.close()
    # Replace original with modified version
    shutil.move(temp_pdf, output_pdf)

    return


def resize_img(source_file,
               width,
               height,
               resolution,
               resize_factor=.92):
    """ Resize an input image or PDF page"""
    # Open image and set dimensions/res
    outer_img = Image(width=width,
                      height=height,
                      resolution=resolution,
                      units='pixelsperinch')

    # if larger than len x wid, fit within box, preserving aspect ratio
    # accordingly determine output dimensions
    tmp_img = Image(filename=source_file)
    tmp_img.transform(resize=f'{width}x{height}>')
    src_wid = int(tmp_img.size[0] * resize_factor)
    src_hgt = int(tmp_img.size[1] * resize_factor)

    del tmp_img

    # Use different libraries to resample PDFs and PNGs
    tmp_source_file = str(source_file)[:-3] + 'temp.png'
    if str(source_file)[-3:] == 'pdf':
        img = Image(filename=source_file)
    else:
        img = PIL_Image.open(source_file)
        # Resize the image using the Lanczos filter (high-quality)
        resized_img = img.resize((src_wid, src_hgt), PIL_Image.LANCZOS)
        # Save the resized image
        resized_img.save(tmp_source_file, dpi=(resolution, resolution))
        del img, resized_img
        img = Image(filename=tmp_source_file)
        Path(tmp_source_file).unlink()

    # set figure positions in page
    left_pos = int((width - img.width) / 2)
    top_pos = int((height - img.height) / 2)
    if '_table_' in Path(source_file).name:
        top_pos = 0
    outer_img.composite(img,
                        left=left_pos,
                        top=top_pos)

    return outer_img


def png_wrapper(png_dir,
                frame,
                output_pdf,
                input_csv,
                page_counter,
                width,
                height,
                resolution,
                top_buffer=0,
                font_color='black'):
    """ Loop through all input images"""
    # Read the CSV data using Pandas
    df = pd.read_csv(input_csv, dtype=str)
    df = df[df['Frame'] == frame]

    # set relevant frame output paths
    frame_outputs = png_dir / frame / 'results'
    frame_gnss_outputs = frame_outputs / f'{frame}_plots'

    # count pages for a given frame ID
    # this includes 3 pages for the standard Calval plots
    frame_pg = 3
    # plus one page for each GNSS vs InSAR comparison
    frame_pg += len([str(i) for i in frame_gnss_outputs.glob('*.png')])
    gnss_input_pngs = [str(i) for i in frame_gnss_outputs.glob('*.png')]

    # loop through each frame page
    for i in range(frame_pg):
        # query VA1 summary map
        if i == 0:
            header_txt = f'VA1 summary for frame {frame}'
            input_pngs = [
                str(i)
                for i in frame_outputs.glob(
                    'Secular_vel_disp_s1_vs_gnss_cartopy_site*.png'
                )
            ][0]

        # query VA1 summary stats
        if i == 1:
            if (
                (df['VA1 (Pass/Fail)'] == 'PASS').iloc[0]
                or (df['VA1 (Pass/Fail)'] == 'FAIL').iloc[0]
            ):
                header_txt = f'VA1 statistics for frame {frame}'
                site_png = [
                    str(i)
                    for i in frame_outputs.glob(
                        'VA1_secular_disp_s1-gnss_velocity_vs_distance_s*.png'
                    )
                ][0]
                tbl_png = [
                    str(i)
                    for i in frame_outputs.glob(
                        'VA1_secular_disp_s1-gnss_velocity_vs_distance_t*.png'
                    )
                ][0]
                input_pngs = site_png + tbl_png
            else:
                continue

        # query VA2 summary stats
        if i == 2:
            if (
                (df['VA2 (Pass/Fail)'] == 'PASS').iloc[0]
                or (df['VA2 (Pass/Fail)'] == 'FAIL').iloc[0]
            ):
                header_txt = f'VA2 statistics for frame {frame}'
                site_png = [
                    str(i)
                    for i in frame_outputs.glob(
                        'VA2_secular_DISP-S1-only_vs_distance_s*.png'
                    )
                ][0]
                tbl_png = [
                    str(i)
                    for i in frame_outputs.glob(
                        'VA2_secular_DISP-S1-only_vs_distance_t*.png'
                    )
                ][0]
                input_pngs = site_png + tbl_png
            else:
                continue

        # query GNSS vs InSAR comparison plots
        if i > 2:
            if (
                (df['VA1 (Pass/Fail)'] == 'PASS').iloc[0]
                or (df['VA1 (Pass/Fail)'] == 'FAIL').iloc[0]
            ):
                header_txt = f'GNSS vs InSAR comparison for frame {frame}'
                gnss_plot_index = i - 3
                input_pngs = [gnss_input_pngs[gnss_plot_index]]
                # skip if there are no GNSS vs InSAR plots
                if input_pngs == []:
                    continue
            else:
                continue

        # build temporary file
        temp_pdf = input_pngs[0].replace('.png', '_temp.pdf')

        # create temp PDF for each image
        for png in input_pngs:
            update_pdf(temp_pdf, png, width, height, resolution)

        if len(input_pngs) > 1:
            # combine images between 2 pages into one
            temp_pdf2 = input_pngs[0].replace('.png', '_fin.pdf')
            merge_pdf_pages(temp_pdf, temp_pdf2,
                            width,
                            height,
                            resolution,
                            top=top_buffer)
        else:
            temp_pdf2 = input_pngs[0].replace('.png', '_temp.pdf')

        # rescale pages
        pdf_rescale = resize_img(temp_pdf2,
                                 width,
                                 height*2,
                                 resolution)
        pdf_rescale.save(filename=temp_pdf2)
        del pdf_rescale

        # add text
        add_text(temp_pdf2, header_txt, font_size=28,
                 font_color=font_color, pos_x=40, pos_y=40)
        # add page number in footer
        page_counter += 1
        ft_x = int((width/2)-10)
        ft_y = int((height*2)-60)
        add_text(temp_pdf2, str(page_counter), font_size=25,
                 bck_color='white', pos_x=ft_x, pos_y=ft_y)

        # append page to cal/val report
        concatenate_pdf_pages(temp_pdf2, output_pdf)

        # delete temp files
        if Path(temp_pdf).exists():
            Path(temp_pdf).unlink()
        if Path(temp_pdf2).exists():
            Path(temp_pdf2).unlink()

    return page_counter


def csv_wrapper(input_csv,
                output_pdf,
                png_dir):
    """ Loop through all input csvs"""
    # Read the CSV data using Pandas
    df = pd.read_csv(input_csv, dtype=str)

    # sort by location and then frame ID
    df = df.sort_values(by=['Frame'])

    # pass frame IDs
    frame_ids = df['Frame'].to_list()
    # get page number for each frame ID
    frame_pg = 2
    pg_num_list = []
    for i in frame_ids:
        # get row corresponding to frame
        df_frame = df[df['Frame'] == i]

        # capture corresponding frame dir for GNSS vs InSAR comparison
        frame_gnss_outputs = png_dir / i / f'results/{i}_plots'

        # capture page number
        pg_num_list.append(frame_pg)

        # count pages for a given frame ID
        # this includes 3 pages for the standard Calval plots
        frame_pg += 3

        # plus one page for each GNSS vs InSAR comparison if VA1 valid
        if (
            (df_frame['VA1 (Pass/Fail)'] == 'PASS').iloc[0]
            or (df_frame['VA1 (Pass/Fail)'] == 'FAIL').iloc[0]
        ):
            frame_pg += len([
                str(i) for i in frame_gnss_outputs.glob('*.png')
            ])
        # remove a page if VA1 not valid
        else:
            frame_pg -= 1

        # remove a page if VA2 not valid
        # plus one page for each GNSS vs InSAR comparison if VA1 valid
        if not (
            (df_frame['VA2 (Pass/Fail)'] == 'PASS').iloc[0]
            or (df_frame['VA2 (Pass/Fail)'] == 'FAIL').iloc[0]
        ):
            frame_pg -= 1

    df['Page #'] = pg_num_list

    # set dims, and scale by expected page size (hardcoded ratios)
    fig_w = 22.2  # 16 * (800/576)
    fig_l = 28.8  # 20.736 * (600/432)

    # Create a figure and axis for the plot
    fig = plt.figure(figsize=(fig_w, fig_l))
    ax = fig.add_subplot(111)

    # Create a table from the Pandas DataFrame
    table = ax.table(cellText=df.values,
                     colLabels=df.columns,
                     loc='center',
                     cellLoc='center')

    # Style the table
    table.auto_set_font_size(False)
    table.set_fontsize(22)
    table.scale(1.2, 2.8)

    # Remove the axis
    ax.axis('off')

    # Save the figure as a PDF
    plt.savefig(output_pdf, pad_inches=0.1)

    # Close the figure
    plt.close()

    # Add table title
    tbl_title = (
        'Summary of CalVal evaluation for all '
        'high priority sites'
    )
    add_text(output_pdf, tbl_title, font_size=28,
             bck_color='white', pos_x=40, pos_y=40)

    # Add caption
    va1_valid_total = (df['VA1 (Pass/Fail)'] == 'PASS').sum()
    va1_valid_total += (df['VA1 (Pass/Fail)'] == 'FAIL').sum()
    passing_rate_va1 = (
        (df['VA1 (Pass/Fail)'] == 'PASS').sum() 
        / va1_valid_total
    ) * 100
    passing_rate_va1 = int(round(passing_rate_va1))
    va2_valid_total = (df['VA2 (Pass/Fail)'] == 'PASS').sum()
    va2_valid_total += (df['VA2 (Pass/Fail)'] == 'FAIL').sum()
    passing_rate_va2 = (
        (df['VA2 (Pass/Fail)'] == 'PASS').sum()
        / va2_valid_total
    ) * 100
    passing_rate_va2 = int(round(passing_rate_va2))
    tbl_txt = f'{passing_rate_va1}% of validation data met VA1 requirement'
    add_text(output_pdf, tbl_txt, font_size=22,
             bck_color='white', pos_x=76, pos_y=1314)
    tbl_txt = f'{passing_rate_va2}% of validation data met VA2 requirement'
    add_text(output_pdf, tbl_txt, font_size=22,
             bck_color='white', pos_x=76, pos_y=1348)
    tbl_txt = (
        '* Frames do not meet temporal sampling requirement (80%) '
        'due to the absence of winter scenes.'
    )
    add_text(output_pdf, tbl_txt, font_size=22,
             bck_color='white', pos_x=76, pos_y=1382)
    tbl_txt = (
        '** Frames do not have sufficient GNSS coverage to assess VA1.'
    )
    add_text(output_pdf, tbl_txt, font_size=22,
             bck_color='white', pos_x=76, pos_y=1416)
    tbl_txt = (
        '*** Frames excluded from VA2 requirement evaluation '
        'as there is significant non-linear motion.'
    )
    add_text(output_pdf, tbl_txt, font_size=22,
             bck_color='white', pos_x=76, pos_y=1450)

    return frame_ids, frame_pg


def update_pdf(output_pdf,
               png,
               width,
               height,
               resolution,
               resize=True,
               resize_factor=.92):
    """ Update a PDF with new content"""
    # initialize output file
    if resize is True:
        pdf_init = resize_img(png,
                              width,
                              height,
                              resolution,
                              resize_factor=resize_factor)
    else:
        pdf_init = Image(filename=png)

    # Save the PDF
    # Initiate output if it does not exist
    if not Path(output_pdf).exists():
        pdf_init.save(filename=output_pdf)
    # otherwise, create a temp file
    else:
        temp_pdf = png.replace('.png', '.pdf')
        # create temp file
        if resize is True:
            pdf_init.save(filename=temp_pdf)
        # re-open cal/val report to append each new page
        concatenate_pdf_pages(temp_pdf, output_pdf)
        # delete temp file
        if resize is True:
            Path(temp_pdf).unlink()

    del pdf_init

    return


def prep_disp_report(inps):
    """ Create a DISP Cal/Val report"""
    # hardcode output page dims/resolution
    width = 1600
    height = 1035
    resolution = 1200

    # pass script home directory
    script_dir = Path(__file__).resolve().parent

    # pass inputs as Path objects
    inps.pdir = Path(inps.pdir)
    inps.outdir = Path(inps.outdir)

    # if necessary create output directory
    if not inps.outdir.exists():
        inps.outdir.mkdir(parents=True, exist_ok=True)

    # hardcode CSV summary name
    csvs_file = script_dir / 'DISP_calval_summary.csv'
    if not csvs_file.exists():
        raise Exception(f'Expected input csv {csvs_file}'
                         'not found in local path')

    # check if specified output file already exists
    output_pdf = inps.outdir / inps.fname
    if output_pdf.exists():
        raise Exception(f'Specified output -f {inps.fname} already exists')

    # initialize output file with validation map
    validation_map = script_dir / 'DISP-S1_Calval_template.pdf'
    output_pdf.write_bytes(Path(validation_map).read_bytes())

    # add table from CSV file
    output_csvpdf = Path(csvs_file).name.replace('.csv', '.pdf')
    output_csvpdf = inps.outdir / output_csvpdf
    # create pdf from table, and pass frame IDs
    frame_ids, frame_pg = csv_wrapper(csvs_file, output_csvpdf,
                                      inps.pdir)
    # add page number in footer
    page_counter = 1
    ft_x = int((width/2)-10)
    ft_y = int((height*2)-50)
    add_text(output_csvpdf, str(page_counter), font_size=25,
             bck_color='white', pos_x=ft_x, pos_y=ft_y)
    # append to report
    if not output_pdf.exists():
        output_pdf.write_bytes(output_csvpdf.read_bytes())
    else:
        concatenate_pdf_pages(output_csvpdf, output_pdf)
    # delete temp file
    output_csvpdf.unlink()

    # if ALE figures captured
    for frame in tqdm(frame_ids, desc="Processing frame"):
        print(frame)
        # loop through each frame png
        page_counter = png_wrapper(inps.pdir, frame, output_pdf, csvs_file,
                                   page_counter, width, height, resolution,
                                   top_buffer=0, font_color='black')

    print(f'PDF created: {output_pdf}')

    return


if __name__ == '__main__':
    main()
    os._exit(0)
