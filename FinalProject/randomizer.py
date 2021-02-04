from random import randint, choice, randrange

id = ['313526675', '321967556', '203192331', '207948381']
questions = []
ranges = [(1, 9), (10, 17), (18, 29), (18, 29), (30, 35), (36, 41)]
for _ in range(6):
    rand_id = choice(id)
    rand_digit = int(rand_id[randint(1, 8)])

    if len(ranges) == 1:
        single_range = ranges.pop(0)
    elif rand_digit % 2 == 0:
        single_range = ranges.pop(randrange(0, len(ranges), 2))
    else:
        single_range = ranges.pop(randrange(1, len(ranges), 2))

    questions.append(randint(single_range[0], single_range[1]))

print(questions)
