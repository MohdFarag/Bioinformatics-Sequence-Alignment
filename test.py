import flet as ft
from utilities import check_sequence

def main(page: ft.Page):
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

    # Sequence input
    def add_sequence_input():
        sequence_input = ft.TextField(label=f"Enter Sequence {len(sequences_inputs_layout.controls)+1}",value="", max_length=40, filled=True)
        sequences_inputs_layout.controls.append(sequence_input)
        page.update()

        return sequence_input

    def num_sequences_minus_click(e):       
        if(int(number_of_sequences_input.value) - 1 < 2):
            return

        number_of_sequences_input.value = str(int(number_of_sequences_input.value) - 1)
        sequences_inputs_layout.controls.pop()
        page.update()
    
    def num_sequences_plus_click(e):
        number_of_sequences_input.value = str(int(number_of_sequences_input.value) + 1)
        add_sequence_input()
        page.update()

    # Number of sequences
    number_of_sequences_div, number_of_sequences_input = create_counter_button(label="Number of sequences", default=2, width=200, height=100, plus_on_click=num_sequences_plus_click, minus_on_click=num_sequences_minus_click)

    # Clear Sequences
    def clear_alignments_action(e):
        sequences_alignments_layout.clean()
        page.update()

    clear_alignments_btn = ft.ElevatedButton("Clear", on_click=clear_alignments_action)

    sequences_inputs_layout = ft.Column([])

    # Default inputs
    for _ in range(20):
        sequence_input = add_sequence_input()

    sequence_text_layout = ft.Column([    
                    number_of_sequences_div,
                    sequences_inputs_layout,
                    ft.Row([
                        clear_alignments_btn
                    ])
                ])

    # Sequences Alignments place
    sequences_alignments_layout = ft.Column(scroll='always',horizontal_alignment= "center")

    primary_page = [
        sequence_text_layout,
        sequences_alignments_layout
    ]

    def pages_routes(route):
        page.views.clear()
        # Primary page
        page.views.append(ft.View('/',primary_page,scroll="always"))


    def view_pop(view):
        page.views.pop()
        my_view = page.views[-1]
        page.go(my_view.route)

    # Settings of the page
    page.title = "Sequence Alignment"
    page.theme_mode = ft.ThemeMode.DARK
    # page.vertical_alignment = ft.MainAxisAlignment.SPACE_EVENLY
    page.scroll = ft.ScrollMode.ALWAYS
    page.auto_scroll = True

    page.on_route_change = pages_routes
    page.on_view_pop = view_pop

    page.go(page.route)
    page.update()


# RUN
ft.app(target=main)
# ft.app(target=main, view=ft.WEB_BROWSER)

