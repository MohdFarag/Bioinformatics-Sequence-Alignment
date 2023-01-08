import flet as ft

from utilities import *
from math import *
from time import sleep

import matplotlib
import matplotlib.pyplot as plt
from flet.matplotlib_chart import MatplotlibChart
matplotlib.use("svg")

color_matrix = []
match_matrix = []
file_path = None
sequences_inputs = []
fig, ax = plt.subplots()

def main(page: ft.Page):
    # Settings of the page
    page.title = "Sequence Alignment"
    page.theme_mode = ft.ThemeMode.DARK
    page.vertical_alignment = ft.MainAxisAlignment.SPACE_EVENLY
    page.scroll = ft.ScrollMode.ALWAYS
    page.auto_scroll = True

    # Counter Button
    def create_counter_button(label, default=0, width=100, height=100, plus_on_click=None, minus_on_click=None):
        def minus_click(e):
            text_field.value = str(int(text_field.value) - 1)
            page.update()

        def plus_click(e):
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

    # matrix plot
    def check_sequence_1(e):
        if not check_sequence(sequence_input_1.value, sequence_type.value):
            sequence_input_1.error_text = f"Invalid {sequence_type.value} Sequence"
            global_alignment_btn.disabled = True
            local_alignment_btn.disabled = True
            msa_btn.disabled = True

            page.update()
            return

        global_alignment_btn.disabled = False
        local_alignment_btn.disabled = False
        msa_btn.disabled = False
        sequence_input_1.error_text = None
        page.update()

    def check_sequence_2(e):
        if not check_sequence(sequence_input_2.value, sequence_type.value):
            sequence_input_2.error_text = f"Invalid {sequence_type.value} Sequence"
            global_alignment_btn.disabled = True
            local_alignment_btn.disabled = True
            msa_btn.disabled = True
            page.update()
            return

        global_alignment_btn.disabled = False
        local_alignment_btn.disabled = False
        msa_btn.disabled = False
        sequence_input_2.error_text = None
        page.update()

    sequence_input_1 = ft.TextField(label="Enter Sequence 1",value="", max_length=40, autofocus=True,filled=True,on_change=check_sequence_1)
    sequence_input_2 = ft.TextField(label="Enter Sequence 2",value="", max_length=40 ,autofocus=True,filled=True,on_change=check_sequence_2)
    
    global sequences_inputs
    sequences_inputs = [sequence_input_1,sequence_input_2]

    sequences_alignments_layout = ft.Column(scroll='always',horizontal_alignment= "center")
    sequence_type = ft.RadioGroup(
        content=ft.Row([
            ft.Radio(value="DNA", label="DNA"),
            ft.Radio(value="RNA", label="RNA"),
            ft.Radio(value="PROTEIN", label="PROTEIN")
        ],
        vertical_alignment="center"))
    sequence_type.value = "DNA" # Default value

    # Scores penalty
    match_div, match_input = create_counter_button("MATCH", 2, 100, 100)
    mismatch_div, mismatch_input = create_counter_button("MISMATCH", -2, 100, 100)
    gap_div, gap_input = create_counter_button("GAP", -1, 100, 100)

    def num_sequences_minus_click(e):
        global sequences_inputs
       
        if int(number_of_sequences_input.value) - 1 > 2:
            msa_btn.visible = True
            global_alignment_btn.visible = False
            local_alignment_btn.visible = False
        else:
            msa_btn.visible = False
            global_alignment_btn.visible = True
            local_alignment_btn.visible = True
        page.update()

        if(int(number_of_sequences_input.value) - 1 < 2):
            return

        number_of_sequences_input.value = str(int(number_of_sequences_input.value) - 1)
        sequences_inputs_layout.controls.remove(sequences_inputs[-1])
        sequences_inputs.pop()
        page.update()
    
    def num_sequences_plus_click(e):
        global sequences_inputs
        number_of_sequences_input.value = str(int(number_of_sequences_input.value) + 1)

        if int(number_of_sequences_input.value) + 1 > 2:
            msa_btn.visible = True
            global_alignment_btn.visible = False
            local_alignment_btn.visible = False
        else:
            msa_btn.visible = False
            global_alignment_btn.visible = True
            local_alignment_btn.visible = True

        def check_sequence_of_input(e):
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

        sequence_input = ft.TextField(label=f"Enter Sequence {len(sequences_inputs)+1}",value="", max_length=40, filled=True,on_change=check_sequence_of_input)
        sequences_inputs.append(sequence_input)
        sequences_inputs_layout.controls.append(sequence_input)
        page.update()

    # Number of sequences
    number_of_sequences_div, number_of_sequences_input = create_counter_button(label="Number of sequences", default=2, width=200, height=100, plus_on_click=num_sequences_plus_click, minus_on_click=num_sequences_minus_click)

    # Score output
    score_output = ft.TextField(label="Score",text_align=ft.TextAlign.CENTER, width=100,disabled=True,filled=True)

    # metrics
    sum_of_pairs_input = ft.TextField(label="sum of pairs",text_align=ft.TextAlign.CENTER, width=100,disabled=True,filled=True)
    mutual_information_input = ft.TextField(label="Mutual information",text_align=ft.TextAlign.CENTER, width=100,disabled=True,filled=True)
    percent_identity_input = ft.TextField(label="percent identity",text_align=ft.TextAlign.CENTER, width=100,disabled=True,filled=True)

    selected_files = ft.Text()
    def pick_files_result(e: ft.FilePickerResultEvent):
        selected_files.value = (", ".join(map(lambda f: f.name, e.files)) if e.files else "Choose file!")

        global file_path
        file_path = (", ".join(map(lambda f: f.path, e.files)) if e.files else None)
        selected_files.update()
    
    pick_files_dialog = ft.FilePicker(on_result=pick_files_result)
    page.overlay.append(pick_files_dialog)
    pick_file_btn = ft.ElevatedButton(
        "Pick files",
        icon=ft.icons.UPLOAD_FILE,
        on_click=lambda _: pick_files_dialog.pick_files(allow_multiple=False)
    )

    # Global Alignment
    def global_alignment_action(e):
        sequences_alignments_layout.controls.clear()
        pb = ft.ProgressBar(width=400)

        sequences_alignments_layout.controls.append(ft.Column([ ft.Text("Process Alignments..."), pb]))
        results = pairwise_global_alignment(sequence_input_1.value, sequence_input_2.value, match_input.value, mismatch_input.value, gap_input.value)
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

        sequences_alignments_layout.controls.clear()
        display_pairwise_alignments(optimal_alignments)

        alignments = []
        for alignment in optimal_alignments:
            for sequence in alignment:
                new_sequence = "".join(sequence)
                alignments.append(new_sequence)

        percent_identity_input.value = percent_identity(alignments)
        sum_of_pairs_input.value = sum_of_pairs(alignments)
        mutual_information_input.value = mutual_information(alignments)
        page.update()
    
    # Local Alignment
    def local_alignment_action(e):
        sequences_alignments_layout.controls.clear()
        pb = ft.ProgressBar(width=400)

        sequences_alignments_layout.controls.append(ft.Column([ ft.Text("Process Alignments..."), pb]))
        results = pairwise_local_alignment(sequence_input_1.value, sequence_input_2.value, match_input.value, mismatch_input.value, gap_input.value)
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

        sequences_alignments_layout.controls.clear()
        display_pairwise_alignments(optimal_alignments)

        alignments = []
        for alignment in optimal_alignments:
            for sequence in alignment:
                new_sequence = "".join(sequence)
                alignments.append(new_sequence)

        percent_identity_input.value = percent_identity(alignments)
        sum_of_pairs_input.value = sum_of_pairs(alignments)
        mutual_information_input.value = mutual_information(alignments)

        page.update() 

    def display_pairwise_alignments(optimal_alignments):
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
           
                sequences_alignments_layout.controls.append(ft.Row(sequence_1))
                sequences_alignments_layout.controls.append(ft.Row(sequence_2))
                sequences_alignments_layout.controls.append(ft.Row([ft.Text(" ")]))

        global fig, ax
        draw_match_matrix(fig,ax,sequence_input_1.value, sequence_input_2.value, match_matrix, color_matrix)

    def display_multiple_alignments(optimal_alignments:dict):
        max_line_size = 40        
        length_of_sequence = len(optimal_alignments[next(iter(optimal_alignments))])
        line_size = min(length_of_sequence,max_line_size)
        last_line_mod = length_of_sequence % line_size
        number_of_lines = ceil((length_of_sequence/max_line_size))

        global fig, ax
        try:
            draw_multiple_sequence_alignment(fig,ax,optimal_alignments)
        except:
            pass
        for k in range(number_of_lines):
            start_index = k*line_size
            end_index = (k+1)*line_size 
            if k == number_of_lines - 1 and k != 0 :
                end_index -= line_size - last_line_mod
            
            for id, alignment in optimal_alignments.items():
                alignment = alignment.decode()
                sequence_line = [ft.Text(fill_word(id,30),width=150)]
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
  
    def select_mode(e):
        if (mode_select.value == "Pick .fasta file"):
            sequence_text_layout.visible = False
            file_layout.visible = True
            plot_btn.visible = True
            show_matrix_btn.visible = False
            msa_btn.visible = True
        else:
            sequence_text_layout.visible = True
            file_layout.visible = False
            plot_btn.visible = False
            show_matrix_btn.visible = True
            if int(number_of_sequences_input.value) > 2:
                msa_btn.visible = True
                global_alignment_btn.visible = False
                local_alignment_btn.visible = False
            else:
                msa_btn.visible = False
                global_alignment_btn.visible = True
                local_alignment_btn.visible = True

        page.update()

    mode_select = ft.Dropdown(
        label="Select option",
        on_change=select_mode,
        width=200,
        height=70,
        options=[
            ft.dropdown.Option("Enter Text Input"),
            ft.dropdown.Option("Pick .fasta file")
    ])

    global_alignment_btn = ft.ElevatedButton("Global Alignment", on_click=global_alignment_action)
    local_alignment_btn = ft.ElevatedButton("Local Alignment", on_click=local_alignment_action)

    # Clear Sequences
    def clear_alignments_action(e):
        sequences_alignments_layout.clean()
        ax.clear()
        score_output.value = None
        page.update()

    def msa_action(e):
        sequences_alignments_layout.controls.clear()
        global file_path
        global sequences_inputs
        if (mode_select.value != "Pick .fasta file"):
            i = 0
            sequences = dict()
            for sequence_input in sequences_inputs:
                sequences[f"Sequence_{i}"] = sequence_input.value
                i += 1

            write_fasta(sequences)
            file_path = output_file
        
        pb = ft.ProgressBar(width=400)

        if file_path != None:
            sequences_alignments_layout.controls.append(ft.Column([ ft.Text("Process Alignments..."), pb]))
            for i in range(0, 101):
                pb.value = i * 0.01
                sleep(0.01)
                page.update()

            results = multiple_sequence_alignment(file_path)
            display_multiple_alignments(results)

            percent_identity_input.value = percent_identity(list(results.values()))
            sum_of_pairs_input.value = sum_of_pairs(list(results.values()))
            mutual_information_input.value = mutual_information(list(results.values()))

        page.update()
        return

    def plot_file():
        draw_multiple_sequence_alignment(fig,ax)
        return
    
    clear_alignments_btn = ft.ElevatedButton("Clear", on_click=clear_alignments_action)
    show_matrix_btn = ft.ElevatedButton("Show Matrix", on_click=lambda _: page.go(f"/plotter/"))

    msa_btn = ft.ElevatedButton("MSA", on_click=msa_action)
    msa_btn.visible = True
    
    plot_btn = ft.ElevatedButton("Plot", on_click=lambda _: page.go(f"/plotter/"))
    plot_btn.visible = False

    grading_layout = ft.Row([
        match_div,
        mismatch_div,
        gap_div
    ])
            
    file_layout = ft.Row([
        pick_file_btn,
        selected_files,
        plot_btn
    ])
    
    analysis_layout=ft.Row([
        score_output,
        sum_of_pairs_input,
        mutual_information_input,
        percent_identity_input
    ])

    global sequences_inputs_layout
    sequences_inputs_layout = ft.Column([sequence_input_1,
                                         sequence_input_2])
    
    sequence_text_layout = ft.Column([    
                    sequence_type,
                    number_of_sequences_div,
                    grading_layout,
                    sequences_inputs_layout,
                    ft.Row([
                        global_alignment_btn,
                        local_alignment_btn,
                        show_matrix_btn,
                        clear_alignments_btn
                    ]) 
                ])

    primary_page = [
        mode_select,
        file_layout,
        sequence_text_layout,
        msa_btn,
        analysis_layout,
        sequences_alignments_layout
    ]

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
        page.views.append(ft.View('/',primary_page))
    
        # matching matrix page
        if page.route == f"/plotter/":
            page.views.append(
                ft.View(
				f"/plotter/",plotter_page)
            )

    def view_pop(view):
        page.views.pop()
        my_view = page.views[-1]
        page.go(my_view.route)

    page.on_route_change = pages_routes
    page.on_view_pop = view_pop

    page.go(page.route)
    page.update()

ft.app(target=main)
# ft.app(target=main, view=ft.WEB_BROWSER)

