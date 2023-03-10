"""
Sequence Alignment - Bioinformatics Project
-----------------------------
    Authors:  Mohamed Ahmed Abdullah & Ismail Tawfik
    Date: 9/1/2023
    Version: 1.0.0
    Status:  Development
    Description:
        GUI for software that can do global & local & multiple sequence alignment.
    Requirements:
        Python 3.9.7
        flet 0.1.0
        matplotlib 3.4.3
        numpy 1.21.2
        bioPython 1.79
"""

from utilities import (
    write_fasta,
    multiple_sequence_alignment,
    pairwise_global_alignment,
    pairwise_local_alignment,
    draw_match_matrix,
    draw_multiple_sequence_alignment,
    check_sequence,
    lists_to_strings,
    fill_word,
    color_of_letter,
    percent_identity,
    mutual_information,
    sum_of_pairs,
    INPUT_LOC)

# Imports
from math import ceil
from time import sleep

# matplotlib imports
import matplotlib
import matplotlib.pyplot as plt
from flet.matplotlib_chart import MatplotlibChart
matplotlib.use("svg")
plt.style.use('dark_background')

# flet imports
import flet as ft

# Global variables
PAIRWISE_LINE_SIZE = 20
MULTIPLE_LINE_SIZE = 40

COLOR_MATRIX = []
MATCH_MATRIX = []
FILE_PATH = None

def main(page: ft.Page):
    # Settings of the page
    page.title = "Sequence Alignment"
    page.theme_mode = ft.ThemeMode.DARK

    fig, axis = plt.subplots()

    # Progress bar download
    def create_progress_bar(label="Processing", width=500, time=101, sleep_time=0.01):
        progress_bar = ft.ProgressBar(width=width)
        sequences_alignments_layout.controls.append(ft.Column(
            [ft.Text(f"{label}..."),
             progress_bar]
        ))
        # Timer
        for i in range(0, time):
            progress_bar.value = i * 0.01
            sleep(sleep_time)
            page.update()

    # Counter Button
    def create_counter_button(label, default=0, width=100, height=100, plus_on_click=None, minus_on_click=None):
        def minus_click(_):
            text_field.value = str(int(text_field.value) - 1)
            page.update()

        def plus_click(_):
            text_field.value = str(int(text_field.value) + 1)
            page.update()

        if plus_on_click is None:
            plus_on_click = plus_click
        
        if minus_on_click is None:
            minus_on_click = minus_click

        minus_button = ft.IconButton(ft.icons.REMOVE, on_click=minus_on_click)
        plus_button = ft.IconButton(ft.icons.ADD, on_click=plus_on_click)
        text_field = ft.TextField(label=str(label), value=str(default), text_align=ft.TextAlign.CENTER, width=width, height=height, filled=True)

        layout = ft.Row([
                        minus_button,
                        text_field,
                        plus_button
                    ], vertical_alignment="center")

        return layout, text_field

    # Sequence input
    def add_sequence_input():
        def check_sequence_of_input(_):
            if not check_sequence(sequence_input.value, sequence_type.value):
                sequence_input.error_text = f"Invalid {sequence_type.value} Sequence"
                global_alignment_btn.disabled = True
                local_alignment_btn.disabled = True
                msa_btn.disabled = True
                page.update()
                return

            global_alignment_btn.disabled = False
            local_alignment_btn.disabled = False
            msa_btn.disabled = False
            sequence_input.error_text = None
            page.update()

        sequence_input = ft.TextField(label=f"Enter Sequence {len(sequences_inputs_layout.controls)+1}",value="", max_length=25, filled=True,on_change=check_sequence_of_input)
        sequences_inputs_layout.controls.append(sequence_input)
        page.update()

        return sequence_input

    # Mode select file input or text input
    def select_mode(_):
        if mode_select.value == "Pick .fasta file":
            page.go("/file/")
        else:
            page.go("/")

        page.update()

    mode_select = ft.Dropdown(
        label="Select option",
        on_change=select_mode,
        width=200,
        height=70,
        value="Enter Text Input",
        options=[
            ft.dropdown.Option("Enter Text Input"),
            ft.dropdown.Option("Pick .fasta file")
        ]
    )

    sequence_type = ft.RadioGroup(
        content=ft.Row([
            ft.Radio(value="DNA", label="DNA"),
            ft.Radio(value="RNA", label="RNA"),
            ft.Radio(value="PROTEIN", label="PROTEIN")
        ],
        vertical_alignment="center"),
        value="DNA")
    
    # Grading penalties
    match_div, match_input = create_counter_button("MATCH", 2, 100, 100)
    mismatch_div, mismatch_input = create_counter_button("MISMATCH", -2, 100, 100)
    gap_div, gap_input = create_counter_button("GAP", -1, 100, 100)

    grading_layout = ft.Row([
        match_div,
        mismatch_div,
        gap_div
    ])

    def num_sequences_minus_click(_):
        if int(number_of_sequences_input.value) - 1 > 2:
            msa_btn.visible = True
            global_alignment_btn.visible = False
            local_alignment_btn.visible = False
        else:
            msa_btn.visible = False
            global_alignment_btn.visible = True
            local_alignment_btn.visible = True
        page.update()

        if int(number_of_sequences_input.value) - 1 < 2:
            return

        number_of_sequences_input.value = str(int(number_of_sequences_input.value) - 1)
        sequences_inputs_layout.controls.pop()
        page.update()
   
    def num_sequences_plus_click(_):
        number_of_sequences_input.value = str(int(number_of_sequences_input.value) + 1)

        if int(number_of_sequences_input.value) + 1 > 2:
            msa_btn.visible = True
            global_alignment_btn.visible = False
            local_alignment_btn.visible = False
        else:
            msa_btn.visible = False
            global_alignment_btn.visible = True
            local_alignment_btn.visible = True

        add_sequence_input()
        page.update()

    # Number of sequences
    number_of_sequences_div, number_of_sequences_input = create_counter_button(label="Number of sequences", default=2, width=200, height=100, plus_on_click=num_sequences_plus_click, minus_on_click=num_sequences_minus_click)

    # Score output
    score_output = ft.TextField(label="Score",text_align=ft.TextAlign.CENTER, width=100,disabled=True,filled=True)

    # Metrics of alignment
    sum_of_pairs_output = ft.TextField(label="Sum Of Pairs",text_align=ft.TextAlign.CENTER, width=100,disabled=True,filled=True)
    mutual_information_output = ft.TextField(label="Mutual Information",text_align=ft.TextAlign.CENTER, width=100,disabled=True,filled=True)
    percent_identity_output = ft.TextField(label="percent Identity",text_align=ft.TextAlign.CENTER, width=100,disabled=True,filled=True)

    # Analysis outputs layout
    analysis_layout = ft.Row([
        score_output,
        sum_of_pairs_output,
        mutual_information_output,
        percent_identity_output
    ])

    # Global Alignment
    def global_alignment_action(_):
        if sequence_input_1.value == "" or sequence_input_2.value == "":
            return
    
        sequences_alignments_layout.controls.clear()

        global MATCH_MATRIX
        global COLOR_MATRIX    
        results = pairwise_global_alignment(sequence_input_1.value, sequence_input_2.value, match_input.value, mismatch_input.value, gap_input.value)
        print(results)

        optimal_alignments = results["alignments"]
        score = results["score"]
        MATCH_MATRIX = results["matrix"]
        COLOR_MATRIX = results["color"]

        create_progress_bar()
        display_pairwise_alignments(optimal_alignments)

        if len(optimal_alignments[0][0]) == 0:
            return

        alignments = lists_to_strings(optimal_alignments[0])
        score_output.value = score
        percent_identity_output.value = percent_identity(alignments)
        sum_of_pairs_output.value = sum_of_pairs(alignments)
        mutual_information_output.value = mutual_information(alignments)
        page.update()
    
    # Local Alignment
    def local_alignment_action(_):
        if sequence_input_1.value == "" or sequence_input_2.value == "":
            return

        sequences_alignments_layout.controls.clear()

        global MATCH_MATRIX
        global COLOR_MATRIX
        results = pairwise_local_alignment(sequence_input_1.value, sequence_input_2.value, match_input.value, mismatch_input.value, gap_input.value)
        print(results)

        optimal_alignments = results["alignments"]
        score = results["score"]
        MATCH_MATRIX = results["matrix"]
        COLOR_MATRIX = results["color"]
        
        create_progress_bar()
        display_pairwise_alignments(optimal_alignments)

        if len(optimal_alignments[0][0]) == 0:
            return

        alignments = lists_to_strings(optimal_alignments[0])
        score_output.value = score
        percent_identity_output.value = percent_identity(alignments)
        sum_of_pairs_output.value = sum_of_pairs(alignments)
        mutual_information_output.value = mutual_information(alignments)

        page.update() 

    def display_pairwise_alignments(optimal_alignments):     
        sequences_alignments_layout.controls.clear()

        for alignment in optimal_alignments:
            result_1 = alignment[0]
            result_2 = alignment[1]

            if len(result_1) == 0 or len(result_2) == 0:
                sequences_alignments_layout.controls.append(ft.Text("No Appropriate alignment found"))
                page.update()
                return
            
            line_size = min(len(result_1), PAIRWISE_LINE_SIZE)
            last_line_mod = len(result_1) % line_size
            number_of_lines = ceil((len(result_1)/PAIRWISE_LINE_SIZE))

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
           
                sequences_alignments_layout.controls.append(ft.Row(sequence_1))
                sequences_alignments_layout.controls.append(ft.Row(sequence_2))
                sequences_alignments_layout.controls.append(ft.Row([ft.Text(" ")]))

        draw_match_matrix(fig, axis, sequence_input_1.value, sequence_input_2.value, MATCH_MATRIX, COLOR_MATRIX)
  
    global_alignment_btn = ft.ElevatedButton("Global Alignment", on_click=global_alignment_action)
    local_alignment_btn = ft.ElevatedButton("Local Alignment", on_click=local_alignment_action)

    # Multiple Sequence Alignment
    def multiple_sequence_alignment_action(_):
        sequences_alignments_layout.controls.clear()

        # Generate a fasta file from the input sequences
        global FILE_PATH
        if (mode_select.value == "Enter Text Input"):
            i = 0
            sequences = dict()
            for sequence_input in sequences_inputs_layout.controls:
                if sequence_input.value == "":
                    return
                else:
                    sequences[f"Sequence-{i}"] = sequence_input.value
                    i += 1

            write_fasta(sequences)
            FILE_PATH = INPUT_LOC
                 
        if FILE_PATH is None:
            return

        results = multiple_sequence_alignment(FILE_PATH)

        create_progress_bar()
        display_multiple_alignments(results)
       
        score_output.value = None
        percent_identity_output.value = percent_identity(list(results.values()))
        sum_of_pairs_output.value = sum_of_pairs(list(results.values()))
        mutual_information_output.value = mutual_information(list(results.values()))

        page.update()

    def display_multiple_alignments(optimal_alignments:dict):
        sequences_alignments_layout.controls.clear()

        length_of_sequence = len(optimal_alignments[next(iter(optimal_alignments))])
        line_size = min(length_of_sequence, MULTIPLE_LINE_SIZE)
        last_line_mod = length_of_sequence % line_size
        number_of_lines = ceil((length_of_sequence/ MULTIPLE_LINE_SIZE))

        # Draw MSA
        draw_multiple_sequence_alignment(fig, axis, optimal_alignments)

        for k in range(number_of_lines):
            start_index = k*line_size
            end_index = (k+1)*line_size
            if k == number_of_lines - 1 and k != 0 :
                end_index -= line_size - last_line_mod
            
            for sequence_id, alignment in optimal_alignments.items():
                alignment = alignment.decode()
                sequence_line = [ft.Text(fill_word(sequence_id,30),width=150)]
                for i in range(start_index,end_index):
                    if alignment[i] == "-":
                        size = 30
                    else:
                        size = 20

                    sequence_line.append(
                        ft.Container(
                                width=30,
                                height=30,
                                bgcolor=color_of_letter(alignment[i]),
                                content=ft.Text(alignment[i],color=ft.colors.WHITE70,size=size),
                                alignment=ft.alignment.center)
                        )
        
                sequences_alignments_layout.controls.append(ft.Row(sequence_line))
            sequences_alignments_layout.controls.append(ft.Text(
                " "
            ))

    msa_btn = ft.ElevatedButton("MSA", on_click=multiple_sequence_alignment_action)

    # Clear Sequences
    def clear_alignments_action(_):
        score_output.value = None
        percent_identity_output.value = None
        sum_of_pairs_output.value = None
        mutual_information_output.value = None
        
        sequences_alignments_layout.clean()
        axis.clear()
        
        page.update()
        
    clear_alignments_btn = ft.ElevatedButton("Clear", on_click=clear_alignments_action)
    show_matrix_btn = ft.ElevatedButton("Show Matrix", on_click=lambda _: page.go(r"/plotter/"))

    sequences_inputs_layout = ft.Column([])
    # Default inputs
    sequence_input_1 = add_sequence_input()
    sequence_input_2 = add_sequence_input()
  
    sequence_text_layout = ft.Column([
                    sequence_type,
                    number_of_sequences_div,
                    grading_layout,
                    sequences_inputs_layout,
                    ft.Row([
                        global_alignment_btn,
                        local_alignment_btn,
                        msa_btn,
                        show_matrix_btn,
                        clear_alignments_btn
                    ])
                ])

    # Sequences Alignments place
    sequences_alignments_layout = ft.Column(scroll='always',horizontal_alignment= "center")

    primary_page = [
        mode_select,
        sequence_text_layout,
        analysis_layout,
        sequences_alignments_layout
    ]

    # Pick file button and file path text field
    selected_files = ft.Text()
    def pick_files_result(e: ft.FilePickerResultEvent):
        selected_files.value = (", ".join(map(lambda f: f.name, e.files)) if e.files else "Choose file!")

        global FILE_PATH
        FILE_PATH = (", ".join(map(lambda f: f.path, e.files)) if e.files else None)
        selected_files.update()    
    pick_files_dialog = ft.FilePicker(on_result=pick_files_result)
    page.overlay.append(pick_files_dialog)
    pick_file_btn = ft.ElevatedButton(
        "Pick files",
        icon=ft.icons.UPLOAD_FILE,
        on_click=lambda _: pick_files_dialog.pick_files(allow_multiple=False),
    )

    plot_btn = ft.ElevatedButton("Plot", on_click=lambda _: page.go(r"/plotter/"))
           
    file_page = [
        mode_select,
        ft.Row([
            pick_file_btn,
            selected_files,
            msa_btn,
            plot_btn,
        ]),
        sequences_alignments_layout
    ]

    # Plotter page
    matching_matrix_plot = MatplotlibChart(fig, expand=True)
    exit_matrix_btn = ft.ElevatedButton("Exit", on_click=lambda _: page.go("/"))
    plotter_page = [
        ft.Container(
            ft.Column([
                matching_matrix_plot,
                exit_matrix_btn
            ]),
        alignment=ft.alignment.center,
        expand=True
    )]

    def pages_routes(route):
        page.views.clear()
        # Primary page
        page.views.append(ft.View('/', primary_page, scroll="adaptive", vertical_alignment = ft.MainAxisAlignment.CENTER))
    
        # Plotter page
        if page.route == r"/plotter/":
            page.views.append(ft.View(page.route, plotter_page))
        
        # Fasta file page
        if page.route == r"/file/":
            page.views.append(ft.View(page.route, file_page, scroll="adaptive"))

    def view_pop(_):
        page.views.pop()
        my_view = page.views[-1]
        page.go(my_view.route)

    page.on_route_change = pages_routes
    page.on_view_pop = view_pop

    page.go(page.route)
    page.update()

# RUN
ft.app(target=main)
# ft.app(target=main, view=ft.WEB_BROWSER)

