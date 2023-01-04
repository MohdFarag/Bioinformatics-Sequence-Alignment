import flet as ft

from utilities import *
from math import *
from time import sleep
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from flet import Page
from flet.matplotlib_chart import MatplotlibChart
from urllib.parse import urlparse
matplotlib.use("svg")

color_matrix = []
match_matrix = []

fig, ax = plt.subplots()
def main(page: ft.Page):
    page.title = "Sequence Alignment"
    youparams = "graph"
    page.theme_mode = ft.ThemeMode.DARK
    page.vertical_alignment = ft.MainAxisAlignment.SPACE_EVENLY
    page.auto_scroll = True

    # matrix plot
    def check_sequence_1(e):
        if not check_sequence(sequence_input_1.value, sequence_type.value):
            sequence_input_1.error_text = f"Invalid {sequence_type.value} Sequence"
            global_alignment_btn.disabled = True
            local_alignment_btn.disabled = True
            page.update()
            return

        global_alignment_btn.disabled = False
        local_alignment_btn.disabled = False
        sequence_input_1.error_text = None
        page.update()

    def check_sequence_2(e):
        if not check_sequence(sequence_input_2.value, sequence_type.value):
            sequence_input_2.error_text = f"Invalid {sequence_type.value} Sequence"
            global_alignment_btn.disabled = True
            local_alignment_btn.disabled = True
            page.update()
            return

        global_alignment_btn.disabled = False
        local_alignment_btn.disabled = False
        sequence_input_2.error_text = None
        page.update()

    sequence_input_1 = ft.TextField(label="Enter Sequence 1",value="GATTACA", autofocus=True,filled=True,on_change=check_sequence_1)
    sequence_input_2 = ft.TextField(label="Enter Sequence 2",value="GTCGACGCA" ,autofocus=True,filled=True,on_change=check_sequence_2)

    Sequences = ft.Column()
    sequence_type = ft.RadioGroup(content=ft.Row([
    ft.Radio(value="DNA", label="DNA"),
    ft.Radio(value="RNA", label="RNA"),
    ft.Radio(value="PROTEIN", label="PROTEIN")]))
    sequence_type.value = "DNA"

    # scores
    match = ft.TextField(label="MATCH",value="2", text_align=ft.TextAlign.CENTER, width=100)
    mismatch = ft.TextField(label="MISMATCH",value="-2", text_align=ft.TextAlign.CENTER, width=120)
    gap = ft.TextField(label="GAP",value="-1", text_align=ft.TextAlign.CENTER, width=100)

    # num of sequences
    num_sequences=ft.TextField(label="Num of Sequences",value="0", text_align=ft.TextAlign.CENTER, width=200)

    # score output
    score_output=ft.TextField(label="score",text_align=ft.TextAlign.CENTER, width=100,disabled=True,filled=True)

    selected_files = ft.Text()

    def pick_files_result(e: ft.FilePickerResultEvent):
        selected_files.value = (
            ", ".join(map(lambda f: f.name, e.files)) if e.files else "Cancelled!"
        )
        selected_files.update()
    pick_files_dialog = ft.FilePicker(on_result=pick_files_result)

    page.overlay.append(pick_files_dialog)

    def global_alignment_action(e):
        Sequences.controls.clear()
        pb = ft.ProgressBar(width=400)

        Sequences.controls.append(ft.Column([ ft.Text("Process Alignments..."), pb]))
        results = pairwise_global_alignment(sequence_input_1.value, sequence_input_2.value, match.value, mismatch.value, gap.value)
        optimal_alignments = results["alignments"]
        score = results["score"]
        score_output.value = score
        global match_matrix
        global color_matrix
        match_matrix = results["matrix"]
        color_matrix = results["color"]
        
        for i in range(0, 101):
            pb.value = i * 0.01
            sleep(0.01)
            page.update()

        Sequences.controls.clear()

        max_line_size = 20        
        for alignment in optimal_alignments:
            result_1 = alignment[0]
            result_2 = alignment[1]
            
            line_size = min(len(result_1),max_line_size)
            last_line_mod = len(result_1) % line_size
            number_of_lines = ceil((len(result_1)/max_line_size))

            for k in range(number_of_lines):
                sequence_1 = []
                sequence_2 = []
                start_index = k*line_size
                end_index = (k+1)*line_size 
                if k == number_of_lines - 1 and k != 0 :
                    end_index -= line_size - last_line_mod
                for i in range(start_index,end_index):
                    if result_1[i] == "-":
                        size = 40
                    else:
                        size = 20

                    sequence_1.append(ft.Container(
                                width=50,
                                height=50,
                                bgcolor=ft.colors.BLACK87,
                                border_radius=100,
                                content=ft.Text(result_1[i],color=ft.colors.WHITE70,size=size),
                                alignment=ft.alignment.center))

                for j in range(start_index,end_index):
                    if result_2[j] == "-":
                        size = 40
                    else:
                        size = 20               
                    
                    sequence_2.append(ft.Container(
                            width=50,
                            height=50,
                            bgcolor=ft.colors.WHITE70,
                            border_radius=100,
                            content=ft.Text(result_2[j],color=ft.colors.BLACK87,size=size),
                            alignment=ft.alignment.center))
           
                Sequences.controls.append(ft.Row(sequence_1))
                Sequences.controls.append(ft.Row(sequence_2))
                Sequences.controls.append(ft.Row())
                Sequences.controls.append(ft.Row())
                Sequences.controls.append(ft.Row())

        global fig, ax
        draw_match_matrix(fig,ax,sequence_input_1.value, sequence_input_2.value,match_matrix,color_matrix)
        page.update()
    
    # Local alignment
    def local_alignment_action(e):
        Sequences.controls.clear()
        pb = ft.ProgressBar(width=400)

        Sequences.controls.append(ft.Column([ ft.Text("Process Alignments..."), pb]))
        results = pairwise_local_alignment(sequence_input_1.value, sequence_input_2.value, match.value, mismatch.value, gap.value)
        optimal_alignments = results["alignments"]
        score = results["score"]
        score_output.value = score

        global match_matrix
        global color_matrix
        match_matrix = results["matrix"]
        color_matrix = results["color"]
        
        for i in range(0, 101):
            pb.value = i * 0.01
            sleep(0.01)
            page.update()

        Sequences.controls.clear()

        max_line_size = 20
        
        for alignment in optimal_alignments:
            result_1 = alignment[0]
            result_2 = alignment[1]
            
            line_size = min(len(result_1),max_line_size)
            last_line_mod = len(result_1) % line_size
            number_of_lines = ceil((len(result_1)/max_line_size))

            for k in range(number_of_lines):
                sequence_1 = []
                sequence_2 = []
                start_index = k*line_size
                end_index = (k+1)*line_size 
                if k == number_of_lines - 1 and k != 0 :
                    end_index -= line_size - last_line_mod
                for i in range(start_index,end_index):
                    if result_1[i] == "-":
                        size = 40
                    else:
                        size = 20

                    sequence_1.append(ft.Container(
                                width=50,
                                height=50,
                                bgcolor=ft.colors.BLACK87,
                                border_radius=100,
                                content=ft.Text(result_1[i],color=ft.colors.WHITE70,size=size),
                                alignment=ft.alignment.center))

                for j in range(start_index,end_index):
                    if result_2[j] == "-":
                        size = 40
                    else:
                        size = 20               
                    
                    sequence_2.append(ft.Container(
                            width=50,
                            height=50,
                            bgcolor=ft.colors.WHITE70,
                            border_radius=100,
                            content=ft.Text(result_2[j],color=ft.colors.BLACK87,size=size),
                            alignment=ft.alignment.center))
           
                Sequences.controls.append(ft.Row(sequence_1))
                Sequences.controls.append(ft.Row(sequence_2))
                Sequences.controls.append(ft.Row())
                Sequences.controls.append(ft.Row())
                Sequences.controls.append(ft.Row())

        global fig, ax
        draw_match_matrix(fig,ax,sequence_input_1.value, sequence_input_2.value,match_matrix,color_matrix)
        page.update() 

    # Clear Sequences
    def clear_alignments_action(e):
        Sequences.clean()
        ax.clear()
        score_output.value=None
        page.update()
   
   # Match count
    def match_minus_click(e):
        match.value = str(int(match.value) - 1)
        page.update()

    def match_plus_click(e):
        match.value = str(int(match.value) + 1)
        page.update()
    
    # MISMATCH COUNT
    def mismatch_minus_click(e):
        mismatch.value = str(int(mismatch.value) - 1)
        page.update()

    def mismatch_plus_click(e):
        mismatch.value = str(int(mismatch.value) + 1)
        page.update()
    
    # GAP COUNT
    def gap_minus_click(e):
        gap.value = str(int(gap.value) - 1)
        page.update()

    def gap_plus_click(e):
        gap.value = str(int(gap.value) + 1)
        page.update()
    
    # sequences count
    def num_sequences_minus_click(e):
        num_sequences.value = str(int(num_sequences.value) - 1)
        num_sequences.update()

    def num_sequences_plus_click(e):
        num_sequences.value = str(int(num_sequences.value) + 1)
        num_sequences.update()
    
    # select sequences
    def select_sequence(e):
        output_text.value = f"selected sequence is:  {sequence_select.value}"
        page.update()
 
    output_text = ft.Text()
    submit_option_btn = ft.ElevatedButton(text="Select", on_click=select_sequence)
    sequence_select = ft.Dropdown(
        label="Choose Sequence",
        width=200,
        height=70,
        options=[
            ft.dropdown.Option("sequence 1")
        ]
    )

    def select_option(e):
        output_text.value = f"selected option is :  {option_select.value}"
        page.update()

    output_text = ft.Text()
    submit_btn = ft.ElevatedButton(text="Select", on_click=select_option)
    option_select = ft.Dropdown(label="Select option",
        width=200,
        height=70,
        options=[
            ft.dropdown.Option("Enter Sequence"),
            ft.dropdown.Option("pick file"),
        ]
    )

    global_alignment_btn = ft.ElevatedButton("Global Alignment", on_click=global_alignment_action)
    local_alignment_btn = ft.ElevatedButton("Local Alignment", on_click=local_alignment_action)
    clear_alignments_btn = ft.ElevatedButton("Clear", on_click=clear_alignments_action)
    show_matrix_btn = ft.ElevatedButton("Show Matrix", on_click=lambda _: page.go(f"/matching_matrix/{youparams}"))
    

    def route_change(route):
        page.views.clear()
        page.views.append(ft.View('/',[
            sequence_type,
            ft.Row(
                [   
                    ft.IconButton(ft.icons.REMOVE, on_click=match_minus_click),
                    match,
                    ft.IconButton(ft.icons.ADD, on_click=match_plus_click),

                    ft.IconButton(ft.icons.REMOVE, on_click=mismatch_minus_click),
                    mismatch,
                    ft.IconButton(ft.icons.ADD, on_click=mismatch_plus_click),

                    ft.IconButton(ft.icons.REMOVE, on_click=gap_minus_click),
                    gap,
                    ft.IconButton(ft.icons.ADD, on_click=gap_plus_click)
                ]
                
            ),ft.Row([option_select,submit_option_btn]),

            ft.Row(
                [
                    ft.ElevatedButton(
                        "Pick files",
                        icon=ft.icons.UPLOAD_FILE,
                        on_click=lambda _: pick_files_dialog.pick_files(
                            allow_multiple=True)),
                            selected_files,
                ]
            ),
            
            ft.Column
            (
                [
                    ft.Row([
                    ft.IconButton(ft.icons.REMOVE, on_click=num_sequences_minus_click),
                    num_sequences,
                    ft.IconButton(ft.icons.ADD, on_click=num_sequences_plus_click)]),
                    sequence_input_1,
                    sequence_input_2,
                    
                    ft.Row
                    (

                    [

                    global_alignment_btn,
                    local_alignment_btn,
                    clear_alignments_btn,

                    ft.Column
                    (
                            [
                                sequence_select
                            ]
                        
                    )
                    ,submit_btn

                    ]
                    ) 
                    ,Sequences,
                    ft.Row([score_output,show_matrix_btn])
                ]),
                    ft.VerticalDivider()
        ]))

        # GET param from home page
        param = page.route
		# THIS IS GET VALUE AFTER /secondpage/THIS RES HERE
        res = urlparse(param).path.split("/")[-1]
        
        matching_matrix_plot = MatplotlibChart(fig,expand=True)
        exit_matrix_btn = ft.ElevatedButton("Exit", on_click=lambda _: page.go("/"))

        if page.route == f"/matching_matrix/{res}":
            page.views.append(
                ft.View(
				f"/matching_matrix/{res}",
				[
                    ft.Container(ft.Column([matching_matrix_plot,exit_matrix_btn]),
                    alignment=ft.alignment.center,
                    expand=True)
                ])
            )

    page.update()
    def view_pop(view):
        page.views.pop()
        myview = page.views[-1]
        page.go(myview.route)

    page.on_route_change = route_change
    page.on_view_pop = view_pop
    page.go(page.route)

ft.app(target=main)
#ft.app(target=main, view=ft.WEB_BROWSER)

